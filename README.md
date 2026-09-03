# piyakcrypt_ctf

**Category:** Crypto

**Flag:** `COMPFEST18{b1as3d_n0nc3_mt_r3c0v3ry_lll_hnp_go_brr_...}`

Recovering an ECDSA private key by exploiting a predictable pseudo-random
generator (Mersenne Twister) combined with a Hidden Number Problem attack
solved via lattice reduction (LLL + Babai).

---

## Table of contents

- [The challenge](#the-challenge)
- [Source code analysis](#source-code-analysis)
- [The vulnerability](#the-vulnerability)
- [Attack plan](#attack-plan)
- [The math](#the-math)
- [Implementation](#implementation)
- [Result](#result)
- [Credits and references](#credits-and-references)

---

## The challenge

A remote service (`nc <ip> <port>`) exposes a menu that allows you to:

```
[1] Show public records      → view the 5 public keys
[2] Show damaged record      → view 40 "gifted" bits of every key
[3] Request a signature      → request an ECDSA signature for a chosen message
[4] Check side room          → no effect
[5] Open data panel          → view 78 "disguised" numbers (max. 9 times)
[6] Submit code               → submit a candidate private key
[7] Exit
```

The server generates 5 ECDSA private keys (**secp256k1** curve, the same one
used by Bitcoin) and allows signing messages with each up to 4 times. The
goal is to recover any one of the 5 private keys and submit it via option 6.

<details>
<summary>See full <code>chall.py</code></summary>

```python
#!/usr/bin/env python3
import os, sys, hashlib, random

P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
A = 0
B = 7
Gx = 55066263022277343669578718895168534326250603453777594175500187360389116729240
Gy = 32670510020758816978083085130507043184471273380659243275938904335757337482424
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

A_HIGH_SIZE = 20
A_LOW_SIZE  = 20
B_SIZE      = 128
UNIT_COUNT  = 5
SIGN_LIMIT  = 4
TABLE_SIZE  = 78
TABLE_LIMIT = 9
SKIP_COUNT  = 622

# ... (see repository for the full file)
```

</details>

---

## Source code analysis

### 1. The curve

```python
P  = 0xFFFF...FEFFFFFC2F
A  = 0
B  = 7
Gx = 55066263022277343669578718895168534326250603453777594175500187360389116729240
Gy = 32670510020758816978083085130507043184471273380659243275938904335757337482424
N  = 0xFFFF...
```

`A=0, B=7` with these specific parameters instantly identify **secp256k1**.

### 2. Private key construction

```python
tag_high = random.getrandbits(A_HIGH_SIZE)   # 20 bits — random.Random (MT19937)
tag_low  = random.getrandbits(A_LOW_SIZE)    # 20 bits — random.Random (MT19937)
piece    = int.from_bytes(os.urandom(D_BYTES), "big")   # 216 bits — OS CSPRNG

secret = (tag_high << (D_SIZE + A_LOW_SIZE)) | (piece << A_LOW_SIZE) | tag_low
pub    = secret * G
```

The 40 bits at the extremes use `random.getrandbits()` — Python's
general-purpose generator, **not cryptographically secure**. The middle
chunk (216 bits) does use `os.urandom()`. Option **2** directly leaks
`tag_high` and `tag_low`.

### 3. The generator leak — option 5

```python
def panel_value(x, pos):
    salt = (0xA5A5A5A5 + pos * 0x6D2B79F5) & MASK32
    bump = (0x9E3779B9 ^ (pos * 0x85EBCA6B)) & MASK32
    y = rol32(x ^ salt, pos * 7 + 3)
    return (y + bump) & MASK32

v = random.getrandbits(32)
print(panel_value(v, pos))
```

`salt` and `bump` **only depend on `pos`**, which the client already knows.
There is no actual secret in the mixing — it is trivially invertible.

### 4. The nonce of each signature

```python
chunk_a = make_piece(random.getrandbits(64), random.getrandbits(64), total_signatures)
chunk_b = int.from_bytes(os.urandom(C_BYTES), "big")
k = (chunk_a << 128) | chunk_b

s = inv_mod(k, N) * (z + r * secret) % N
```

The top 128 bits of the nonce `k` come once again from the same insecure
generator (`random`); the bottom 128 bits are secure (`os.urandom`).

---

## The vulnerability

A single design flaw shows up in two places in the program:

> **`random.getrandbits()` uses Mersenne Twister (MT19937) — a
> deterministic, publicly invertible algorithm, unsuitable for
> cryptographic use.**

If **624 consecutive 32-bit outputs** of the generator are observed, its
entire internal state (624 × 32 = 19,937 bits) can be reconstructed, and
any future output can be predicted with absolute certainty.

Option 5 is the only window the program leaves open into that shared
generator — and although it tries to disguise the numbers, the disguise is
reversible.

---

## Attack plan

```
1. Request option 5 nine times              → 702 disguised numbers
2. Undo panel_value()                        → 702 raw random outputs
3. Reconstruct MT19937 state                 → from the first 624
4. Predict chunk_a for each future signature → before requesting it
5. Request 4 signatures for each of the 5 units
6. For each unit:
     - set up the HNP system (128 known bits of k, 128 unknown)
     - eliminate the private key by combining pairs of signatures
     - build the lattice basis
     - reduce with LLL
     - find the closest vector with Babai's algorithm
     - derive the candidate private key
7. Verify against the gifted bits (option 2) and the public key (option 1)
8. Submit the key via option 6 → flag
```

---

## The math

### Splitting the nonce into known + unknown

$$k = A + x, \qquad A \text{ known (predicted)}, \quad 0 \le x < 2^{128}$$

### Rewriting the signature equation

$$s \equiv k^{-1}(z + rd) \pmod N \ \Longrightarrow\ x \equiv \underbrace{(s^{-1}z - A)}_{c} + \underbrace{(s^{-1}r)}_{t}\cdot d \pmod N$$

### Eliminating `d` by combining two signatures

$$d \equiv t_0^{-1}(x_0-c_0) \ \Longrightarrow\ x_i - \underbrace{(t_i t_0^{-1})}_{w_i}\, x_0 \equiv \underbrace{c_i - w_i c_0}_{b_i} \pmod N$$

### The lattice (Hidden Number Problem, Boneh–Venkatesan)

$$
B =
\begin{pmatrix}
N & 0 & 0 & 0 \\
0 & N & 0 & 0 \\
0 & 0 & N & 0 \\
w_1 & w_2 & w_3 & 1
\end{pmatrix}
\qquad
\text{target} = (b_1, b_2, b_3, 0)
$$

LLL reduces the basis; Babai's algorithm finds the lattice vector closest
to the target. Its last coordinate is exactly `x₀` — the missing piece of
the nonce.

### Recovering and verifying the key

$$d \equiv t_0^{-1}(x_0-c_0) \pmod N \qquad\qquad d\cdot G \stackrel{?}{=} \text{pub}$$

---

## Implementation

Pure Python solver, no `pwntools` (avoids heavy build dependencies like
`unicorn`) — only the standard library's `socket`.

**Dependencies:**
```bash
pip install randcrack fpylll cysignals
```

**Key pieces of the solver:**

| Function | Role |
|---|---|
| `desdisfraza_valor()` | Inverts `panel_value()` |
| `clonar_random()` | Reconstructs the MT19937 state with `randcrack` |
| `predecir_chunk_a()` | Predicts the top 128 bits of the nonce |
| `buscar_clave_unidad()` | Sets up the HNP, builds the lattice, runs LLL + Babai |
| `hacer_lll()` | LLL reduction via `fpylll` |
| `babai_vecino_mas_cercano()` | CVP via Babai's algorithm (exact arithmetic with `fractions.Fraction`) |

The full code is in [`solve.py`](./solve.py).

---

## Result

```
$ sage -python solve.py
[*] Collecting data-panel leaks (9x)...
    ...
[*] Cracking MT19937 state...
    alignment check: OK
[*] Reading damaged record (sanity-check bits of every secret)...
    tag_high=0x87fdf tag_low=0x081f7
[*] Reading public keys...
    got 5 public keys
[*] Requesting signatures for all units...
    ...
[*] Running HNP lattice attack per unit...
[+] Recovered secret for unit 0: 61510954339560518939547126329213547641079480598258146797275258950427988754935
[*] Submitting...

  Accepted for unit #0.
```

### 🚩 Flag

```
COMPFEST18{b1as3d_n0nc3_mt_r3c0v3ry_lll_hnp_go_brr_727e3a9724b244c1}
```

---

## Credits and references

The underlying technique (partially leaked ECDSA nonce + Hidden Number
Problem solved via lattice reduction) has been public since 1996 and has
been reused in numerous CTF challenges. I recognized the pattern while
comparing with:

- **Sign Wars — SECCON CTF 2021** ([writeup](https://org.anize.rs/SECCON-2021/crypto/signwars)) —
  a challenge with the same family of vulnerability (crackable MT19937 +
  biased ECDSA nonce). Their setup had one additional unknown (the signed
  message was also hidden), so they used Z3 plus an external CVP solver;
  in this challenge the message is known, so the system could be solved
  with a lighter, self-written implementation of Babai/LLL.
- **D. Boneh, R. Venkatesan**, *"Hardness of Computing the Most Significant
  Bits of Secret Keys in Diffie-Hellman and Related Schemes,"* CRYPTO 1996 —
  the academic origin of the Hidden Number Problem.
- **J. Breitner, N. Heninger**, *"Biased Nonce Sense: Lattice Attacks
  against Weak ECDSA Signatures in Cryptocurrencies,"* 2019
  ([eprint](https://eprint.iacr.org/2019/023.pdf)) — formalization of the
  lattice attack on partially known ECDSA nonces.
- **"Cracking Random Number Generators"** series by jazzy.id.au — the
  classic reference on Mersenne Twister state recovery.
- [`randcrack`](https://github.com/tna0y/Python-random-module-cracker) and
  [`fpylll`](https://github.com/fplll/fpylll) — libraries used in the
  implementation.


---
