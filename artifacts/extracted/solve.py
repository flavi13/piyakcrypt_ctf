#!/usr/bin/env python3
import re
import socket
import sys
from fractions import Fraction as Fr

from randcrack import RandCrack
from fpylll import IntegerMatrix, LLL

HOST = "34.2.147.230"
PORT = 3002


P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
A = 0
B = 7
Gx = 55066263022277343669578718895168534326250603453777594175500187360389116729240
Gy = 32670510020758816978083085130507043184471273380659243275938904335757337482424
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

A_HIGH_SIZE = 20
A_LOW_SIZE = 20
B_SIZE = 128
UNIT_COUNT = 5
SIGN_LIMIT = 4
TABLE_SIZE = 78
TABLE_LIMIT = 9
MASK32 = (1 << 32) - 1
MASK64 = (1 << 64) - 1
D_SIZE = 256 - A_HIGH_SIZE - A_LOW_SIZE


def inv_mod(x, m):
    return pow(x, -1, m)


def ec_add(P1, P2):
    if P1 is None:
        return P2
    if P2 is None:
        return P1
    x1, y1 = P1
    x2, y2 = P2
    if x1 == x2 and (y1 + y2) % P == 0:
        return None
    if P1 == P2:
        lam = (3 * x1 * x1 + A) * inv_mod(2 * y1 % P, P) % P
    else:
        lam = (y2 - y1) * inv_mod((x2 - x1) % P, P) % P
    x3 = (lam * lam - x1 - x2) % P
    y3 = (lam * (x1 - x3) - y1) % P
    return (x3, y3)


def ec_mul(k, pt):
    if k % N == 0 or pt is None:
        return None
    k %= N
    R = None
    Q = pt
    while k:
        if k & 1:
            R = ec_add(R, Q)
        Q = ec_add(Q, Q)
        k >>= 1
    return R


def rol32(x, r):
    r &= 31
    return ((x << r) | (x >> (32 - r))) & MASK32


def ror32(x, r):
    r &= 31
    return ((x >> r) | (x << (32 - r))) & MASK32


def rol64(x, r):
    r &= 63
    return ((x << r) | (x >> (64 - r))) & MASK64


def panel_value(x, pos):
    salt = (0xA5A5A5A5 + pos * 0x6D2B79F5) & MASK32
    bump = (0x9E3779B9 ^ (pos * 0x85EBCA6B)) & MASK32
    y = rol32(x ^ salt, pos * 7 + 3)
    return (y + bump) & MASK32


def unpanel_value(v, pos):
    salt = (0xA5A5A5A5 + pos * 0x6D2B79F5) & MASK32
    bump = (0x9E3779B9 ^ (pos * 0x85EBCA6B)) & MASK32
    y = (v - bump) & MASK32
    x = ror32(y, pos * 7 + 3)
    return x ^ salt


def fold_piece(x, pos, lane):
    x ^= ((pos + 1) * 0xD6E8FEB86659FD93 + lane * 0xA0761D6478BD642F) & MASK64
    x = rol64(x, 17 + pos * 9 + lane * 23)
    x = (x * 0x9E6C63D0676A9A99 + 0xD1B54A32D192ED03) & MASK64
    return x


def make_piece(a, b, pos):
    return (fold_piece(a, pos, 0) << 64) | fold_piece(b, pos, 1)


def crack_rng(panel_dict):
    n = len(panel_dict)
    words = [unpanel_value(panel_dict[p], p) for p in range(n)]
    rc = RandCrack()
    for w in words[:624]:
        rc.submit(w)
    ok = True
    for p in range(624, n):
        if rc.predict_getrandbits(32) != words[p]:
            ok = False
    return rc, ok


def predict_chunk_a(rc, pos):
    a = rc.predict_getrandbits(64) # flavii: El predict no es seguro criptográficamente
    b = rc.predict_getrandbits(64)
    return make_piece(a, b, pos)


def gram_schmidt(Bmat):
    n = len(Bmat)
    dim = len(Bmat[0])
    Bstar = [None] * n
    mu = [[Fr(0)] * n for _ in range(n)]
    for i in range(n):
        Bstar[i] = [Fr(int(x)) for x in Bmat[i]]
        for j in range(i):
            num = sum(Fr(int(Bmat[i][k])) * Bstar[j][k] for k in range(dim))
            den = sum(Bstar[j][k] * Bstar[j][k] for k in range(dim))
            mu[i][j] = num / den
            Bstar[i] = [Bstar[i][k] - mu[i][j] * Bstar[j][k] for k in range(dim)]
    return Bstar


def babai_nearest_plane(Bmat, target):
    n = len(Bmat)
    dim = len(target)
    Bstar = gram_schmidt(Bmat)
    b = [Fr(int(x)) for x in target]
    coeffs = [0] * n
    for i in reversed(range(n)):
        num = sum(b[k] * Bstar[i][k] for k in range(dim))
        den = sum(Bstar[i][k] * Bstar[i][k] for k in range(dim))
        c = round(num / den) if den != 0 else 0
        coeffs[i] = c
        b = [b[k] - c * Fr(int(Bmat[i][k])) for k in range(dim)]
    closest = [0] * dim
    for i in range(n):
        for k in range(dim):
            closest[k] += coeffs[i] * int(Bmat[i][k])
    return closest


def lll_reduce(Bmat):
    n = len(Bmat)
    dim = len(Bmat[0])
    M = IntegerMatrix(n, dim)
    for i in range(n):
        for j in range(dim):
            M[i, j] = int(Bmat[i][j])
    LLL.reduction(M)

    return [[int(M[i, j]) for j in range(dim)] for i in range(n)]


def solve_unit_secret(sigs):
    """
    sigs: list of (r, s, z, A) for signatures made with ONE secret key d,
          A = (predicted_top_128_bits_of_k << 128) % N.
    Yields candidate d values to check against the known public key.
    """
    m = len(sigs)
    if m < 2:
        return

    t, c = [], []
    for (r, s, z, Ai) in sigs:
        sinv = pow(s, -1, N)
        t.append((sinv * r) % N)
        c.append((sinv * z - Ai) % N)

    t0inv = pow(t[0], -1, N)
    w, b = [], []
    for i in range(1, m):
        wi = (t[i] * t0inv) % N
        bi = (c[i] - wi * c[0]) % N
        w.append(wi)
        b.append(bi)

    k = len(w)
    dim = k + 1
    Bmat = []
    for i in range(k):
        row = [0] * dim
        row[i] = N
        Bmat.append(row)
    Bmat.append(w + [1])

    Bred = lll_reduce(Bmat)
    target = b + [0]
    closest = babai_nearest_plane(Bred, target)
    base_alpha = closest[-1]

    candidates = set()
    for sign in (1, -1):
        for delta in range(-4, 5):
            candidates.add(sign * base_alpha + delta)

    for alpha in candidates:
        d = (t0inv * (alpha - c[0])) % N
        yield d

class Conn:
    def __init__(self, host, port, timeout=10):
        self.sock = socket.create_connection((host, port), timeout=timeout)
        self.buf = b""

    def _fill(self):
        chunk = self.sock.recv(65536)
        if not chunk:
            raise ConnectionError("remote closed the connection")
        self.buf += chunk

    def recvuntil(self, marker, decode=True):
        while marker not in self.buf:
            self._fill()
        idx = self.buf.index(marker) + len(marker)
        data, self.buf = self.buf[:idx], self.buf[idx:]
        return data.decode(errors="replace") if decode else data

    def sendline(self, data):
        if isinstance(data, str):
            data = data.encode()
        self.sock.sendall(data + b"\n")

    def recvall(self, timeout=5):
        self.sock.settimeout(timeout)
        try:
            while True:
                self._fill()
        except (socket.timeout, ConnectionError):
            pass
        out, self.buf = self.buf, b""
        return out.decode(errors="replace")

    def close(self):
        try:
            self.sock.close()
        except OSError:
            pass


def get_panel(io):
    io.sendline(b"5")
    blob = io.recvuntil(b"menu> ")
    out = {}
    for mobj in re.finditer(r"entry_(\d+)\s*=\s*0x([0-9a-fA-F]+)", blob):
        pos = int(mobj.group(1))
        val = int(mobj.group(2), 16)
        out[pos] = val
    return out


def get_damaged_record(io):
    io.sendline(b"2")
    blob = io.recvuntil(b"menu> ")
    th = int(re.search(r"record_high\s*=\s*0x([0-9a-fA-F]+)", blob).group(1), 16)
    tl = int(re.search(r"record_low\s*=\s*0x([0-9a-fA-F]+)", blob).group(1), 16)
    return th, tl


def get_public_keys(io):
    io.sendline(b"1")
    blob = io.recvuntil(b"menu> ")
    pubs = {}
    for mobj in re.finditer(
        r"Unit #(\d+):\s*X = 0x([0-9a-fA-F]+)\s*Y = 0x([0-9a-fA-F]+)", blob
    ):
        idx = int(mobj.group(1))
        x = int(mobj.group(2), 16)
        y = int(mobj.group(3), 16)
        pubs[idx] = (x, y)
    return pubs


def request_signature(io, unit_idx, msg):
    io.sendline(b"3")
    io.recvuntil(b"Choose unit")
    io.sendline(str(unit_idx).encode())
    io.recvuntil(b"Message")
    io.sendline(msg)
    blob = io.recvuntil(b"menu> ")
    z = int(re.search(r"z\s*=\s*(\d+)", blob).group(1))
    r = int(re.search(r"r\s*=\s*(\d+)", blob).group(1))
    s = int(re.search(r"s\s*=\s*(\d+)", blob).group(1))
    return z, r, s


def submit_code(io, guess):
    io.sendline(b"6")
    io.recvuntil(b"Code")
    io.sendline(str(guess).encode())
    return io.recvall(timeout=5)


def main():
    io = Conn(HOST, PORT)
    io.recvuntil(b"menu> ")  

    print("[*] Collecting data-panel leaks (9x)...")
    panel = {}
    for i in range(TABLE_LIMIT):
        panel.update(get_panel(io))
        print(f"    read batch {i+1}/{TABLE_LIMIT}, total entries: {len(panel)}")

    print("[*] Cracking MT19937 state...")
    rc, aligned = crack_rng(panel)
    print("    alignment check:", "OK" if aligned else "MISMATCH (attack may still work if only tail entries mismatch)")

    print("[*] Reading damaged record (sanity-check bits of every secret)...")
    tag_high, tag_low = get_damaged_record(io)
    print(f"    tag_high=0x{tag_high:05x} tag_low=0x{tag_low:05x}")

    print("[*] Reading public keys...")
    pubs = get_public_keys(io)
    print(f"    got {len(pubs)} public keys")

    print("[*] Requesting signatures for all units...")
    sig_data = {i: [] for i in range(UNIT_COUNT)}
    total_sigs = 0
    for unit_idx in range(UNIT_COUNT):
        for j in range(SIGN_LIMIT):
            pos = total_sigs  
            chunk_a_pred = predict_chunk_a(rc, pos)
            msg = f"solve-{total_sigs}".encode()
            z, r, s = request_signature(io, unit_idx, msg)
            A_i = (chunk_a_pred << (256 - B_SIZE)) % N
            sig_data[unit_idx].append((r, s, z, A_i))
            total_sigs += 1
            print(f"    unit {unit_idx} sig {j+1}/{SIGN_LIMIT} (total {total_sigs})")

    print("[*] Running HNP lattice attack per unit...")
    found_secret = None
    found_unit = None
    for unit_idx, sigs in sig_data.items():
        pub = pubs.get(unit_idx)
        if pub is None:
            continue
        for d in solve_unit_secret(sigs):
            if (d >> (D_SIZE + A_LOW_SIZE)) != tag_high:
                continue
            if (d & ((1 << A_LOW_SIZE) - 1)) != tag_low:
                continue
            if ec_mul(d, (Gx, Gy)) == pub:
                found_secret = d
                found_unit = unit_idx
                break
        if found_secret is not None:
            break

    if found_secret is None:
        print("[-] Failed to recover any unit's secret. Dumping debug info below.")
        print("    sig_data:", sig_data)
        io.close()
        sys.exit(1)

    print(f"[+] Recovered secret for unit {found_unit}: {found_secret}")
    print("[*] Submitting...")
    result = submit_code(io, found_secret)
    print(result)


if __name__ == "__main__":
    main()
