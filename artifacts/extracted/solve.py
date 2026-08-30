#!/usr/bin/env python3

import re
import socket
import sys #parar programa en caso de error
from fractions import Fraction as Fr #resultado exacto sin floats

from randcrack import RandCrack # clonar el random de python (algoritmo para saber que saldra despues)
from fpylll import IntegerMatrix, LLL # hacer la reduccion 

IP = "34.2.147.230"
PUERTO = 3002

# los saque del codigo fuente 
P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
A = 0
B = 7
Gx = 55066263022277343669578718895168534326250603453777594175500187360389116729240
Gy = 32670510020758816978083085130507043184471273380659243275938904335757337482424
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

# estas constantes son del propio reto, copiadas del .py que nos dieron
A_HIGH_SIZE = 20
A_LOW_SIZE = 20
B_SIZE = 128
UNIT_COUNT = 5
SIGN_LIMIT = 4       # firmas maximas por unidad
TABLE_SIZE = 78
TABLE_LIMIT = 9       # veces que se puede pedir la opcion 5
MASK32 = (1 << 32) - 1
MASK64 = (1 << 64) - 1
D_SIZE = 256 - A_HIGH_SIZE - A_LOW_SIZE


def inv_mod(x, m):

    return pow(x, -1, m)


# suma de puntos
def suma_puntos(p1, p2):
    if p1 is None:
        return p2
    if p2 is None:
        return p1
    x1, y1 = p1
    x2, y2 = p2
    if x1 == x2 and (y1 + y2) % P == 0:
        return None
    if p1 == p2:
        lam = (3 * x1 * x1 + A) * inv_mod(2 * y1 % P, P) % P
    else:
        lam = (y2 - y1) * inv_mod((x2 - x1) % P, P) % P
    x3 = (lam * lam - x1 - x2) % P
    y3 = (lam * (x1 - x3) - y1) % P
    return (x3, y3)


def mult_punto(k, pt):
    # multiplicacion 
    if k % N == 0 or pt is None:
        return None
    k %= N
    res = None
    q = pt
    while k:
        if k & 1:
            res = suma_puntos(res, q)
        q = suma_puntos(q, q)
        k >>= 1
    return res

# reconstruir el pedazo de nonce que predecimos con randcrack

def rol32(x, r):
    r &= 31
    return ((x << r) | (x >> (32 - r))) & MASK32


def ror32(x, r):
    # rotar al reves, para deshacer el rol32 de arriba
    r &= 31
    return ((x >> r) | (x << (32 - r))) & MASK32


def rol64(x, r):
    r &= 63
    return ((x << r) | (x >> (64 - r))) & MASK64


def disfraza_valor(x, pos):
    #  panel_value() del reto 
    salt = (0xA5A5A5A5 + pos * 0x6D2B79F5) & MASK32
    bump = (0x9E3779B9 ^ (pos * 0x85EBCA6B)) & MASK32
    y = rol32(x ^ salt, pos * 7 + 3)
    return (y + bump) & MASK32


def desdisfraza_valor(v, pos):
    # funcion inversa
    salt = (0xA5A5A5A5 + pos * 0x6D2B79F5) & MASK32
    bump = (0x9E3779B9 ^ (pos * 0x85EBCA6B)) & MASK32
    y = (v - bump) & MASK32
    x = ror32(y, pos * 7 + 3)
    return x ^ salt


def revuelve_pieza(x, pos, lane): 
    x ^= ((pos + 1) * 0xD6E8FEB86659FD93 + lane * 0xA0761D6478BD642F) & MASK64
    x = rol64(x, 17 + pos * 9 + lane * 23)
    x = (x * 0x9E6C63D0676A9A99 + 0xD1B54A32D192ED03) & MASK64
    return x


def arma_trozo_nonce(a, b, pos):
    # make_piece() del reto, junta dos numeros de 64 bits en uno de 128
    return (revuelve_pieza(a, pos, 0) << 64) | revuelve_pieza(b, pos, 1)

#

def clonar_random(panel):
    n = len(panel)
    numeros = [desdisfraza_valor(panel[p], p) for p in range(n)]

    rc = RandCrack()
    for w in numeros[:624]:
        rc.submit(w)

    #  comprobación 
    todo_bien = True
    for p in range(624, n):
        if rc.predict_getrandbits(32) != numeros[p]:
            todo_bien = False

    return rc, todo_bien


def predecir_chunk_a(rc, pos):
    # cada firma gasta 2 llamadas a getrandbits(64) para construir
    a = rc.predict_getrandbits(64)
    b = rc.predict_getrandbits(64)
    return arma_trozo_nonce(a, b, pos)



# LLL

def gram_schmidt(matriz):
    n = len(matriz)
    dim = len(matriz[0])
    base_ortog = [None] * n
    mu = [[Fr(0)] * n for _ in range(n)]
    for i in range(n):
        base_ortog[i] = [Fr(int(x)) for x in matriz[i]]
        for j in range(i):
            num = sum(Fr(int(matriz[i][k])) * base_ortog[j][k] for k in range(dim))
            den = sum(base_ortog[j][k] * base_ortog[j][k] for k in range(dim))
            mu[i][j] = num / den
            base_ortog[i] = [base_ortog[i][k] - mu[i][j] * base_ortog[j][k] for k in range(dim)]
    return base_ortog


def babai_vecino_mas_cercano(matriz, objetivo):
    n = len(matriz)
    dim = len(objetivo)
    base_ortog = gram_schmidt(matriz)
    b = [Fr(int(x)) for x in objetivo]
    coefs = [0] * n
    for i in reversed(range(n)):
        num = sum(b[k] * base_ortog[i][k] for k in range(dim))
        den = sum(base_ortog[i][k] * base_ortog[i][k] for k in range(dim))
        c = round(num / den) if den != 0 else 0
        coefs[i] = c
        b = [b[k] - c * Fr(int(matriz[i][k])) for k in range(dim)]

    cercano = [0] * dim
    for i in range(n):
        for k in range(dim):
            cercano[k] += coefs[i] * int(matriz[i][k])
    return cercano


def hacer_lll(matriz):
    #  libreria fpylll, que ya tiene el LLL implementado
    n = len(matriz)
    dim = len(matriz[0])
    m = IntegerMatrix(n, dim)
    for i in range(n):
        for j in range(dim):
            m[i, j] = int(matriz[i][j])
    LLL.reduction(m)
    #  forzar int() por problema con sagemath en mi consola
    return [[int(m[i, j]) for j in range(dim)] for i in range(n)]


def buscar_clave_unidad(firmas):
   
    m = len(firmas)
    if m < 2:
        # con 1 sola firma no se puede montar el sistema, se necesitan 2
        return

    t_list, c_list = [], []
    for (r, s, z, Ai) in firmas:
        s_inv = pow(s, -1, N)
        t_list.append((s_inv * r) % N)
        c_list.append((s_inv * z - Ai) % N)

    # usamos la firma 0 como "base" para cancelar la clave d de las demas
    t0_inv = pow(t_list[0], -1, N)
    w_list, b_list = [], []
    for i in range(1, m):
        wi = (t_list[i] * t0_inv) % N
        bi = (c_list[i] - wi * c_list[0]) % N
        w_list.append(wi)
        b_list.append(bi)

    dim_extra = len(w_list)
    dim = dim_extra + 1

    # construimos la matriz para el LLL
    matriz = []
    for i in range(dim_extra):
        fila = [0] * dim
        fila[i] = N
        matriz.append(fila)
    matriz.append(w_list + [1])

    matriz_reducida = hacer_lll(matriz)
    objetivo = b_list + [0]
    cercano = babai_vecino_mas_cercano(matriz_reducida, objetivo)

    x0_candidato = cercano[-1]
    candidatos = set()
    for signo in (1, -1):
        for delta in range(-4, 5):
            candidatos.add(signo * x0_candidato + delta)

    for x0 in candidatos:
        d = (t0_inv * (x0 - c_list[0])) % N
        yield d

class Conexion:
    def __init__(self, host, puerto, timeout=10):
        self.sock = socket.create_connection((host, puerto), timeout=timeout)
        self.buffer = b""

    def _leer_mas(self):
        trozo = self.sock.recv(65536)
        if not trozo:
            raise ConnectionError("el server ha cerrado la conexion")
        self.buffer += trozo

    def recvuntil(self, marca):
        # va leyendo hasta que encuentra el texto que buscamos
        while marca not in self.buffer:
            self._leer_mas()
        idx = self.buffer.index(marca) + len(marca)
        data, self.buffer = self.buffer[:idx], self.buffer[idx:]
        return data.decode(errors="replace")

    def sendline(self, data):
        if isinstance(data, str):
            data = data.encode()
        self.sock.sendall(data + b"\n")

    def recvall(self, timeout=5):
        self.sock.settimeout(timeout)
        try:
            while True:
                self._leer_mas()
        except (socket.timeout, ConnectionError):
            pass
        out, self.buffer = self.buffer, b""
        return out.decode(errors="replace")

    def close(self):
        try:
            self.sock.close()
        except OSError:
            pass


# funciones para el Menu

def pedir_panel(io):
    io.sendline(b"5")
    texto = io.recvuntil(b"menu> ")
    salida = {}
    for m in re.finditer(r"entry_(\d+)\s*=\s*0x([0-9a-fA-F]+)", texto):
        pos = int(m.group(1))
        val = int(m.group(2), 16)
        salida[pos] = val
    return salida


def pedir_record_dañado(io):
    io.sendline(b"2")
    texto = io.recvuntil(b"menu> ")
    th = int(re.search(r"record_high\s*=\s*0x([0-9a-fA-F]+)", texto).group(1), 16)
    tl = int(re.search(r"record_low\s*=\s*0x([0-9a-fA-F]+)", texto).group(1), 16)
    return th, tl


def pedir_claves_publicas(io):
    io.sendline(b"1")
    texto = io.recvuntil(b"menu> ")
    pubs = {}
    for m in re.finditer(r"Unit #(\d+):\s*X = 0x([0-9a-fA-F]+)\s*Y = 0x([0-9a-fA-F]+)", texto):
        idx = int(m.group(1))
        x = int(m.group(2), 16)
        y = int(m.group(3), 16)
        pubs[idx] = (x, y)
    return pubs


def pedir_firma(io, unidad, mensaje):
    io.sendline(b"3")
    io.recvuntil(b"Choose unit")
    io.sendline(str(unidad).encode())
    io.recvuntil(b"Message")
    io.sendline(mensaje)
    texto = io.recvuntil(b"menu> ")
    z = int(re.search(r"z\s*=\s*(\d+)", texto).group(1))
    r = int(re.search(r"r\s*=\s*(\d+)", texto).group(1))
    s = int(re.search(r"s\s*=\s*(\d+)", texto).group(1))
    return z, r, s


def mandar_codigo(io, clave):
    io.sendline(b"6")
    io.recvuntil(b"Code")
    io.sendline(str(clave).encode())
    return io.recvall(timeout=5)


def main():
    io = Conexion(IP, PUERTO)
    io.recvuntil(b"menu> ")  
    print("[*] pidiendo el panel de datos 9 veces...")
    panel = {}
    for i in range(TABLE_LIMIT):
        panel.update(pedir_panel(io))
        print(f"    van {i+1}/{TABLE_LIMIT}, total numeros: {len(panel)}")

    print("[*] clonando el generador random...")
    rc, ok = clonar_random(panel)
    print("    check:", "todo cuadra" if ok else "algo no cuadra pero seguimos igualmente")

    print("[*] leyendo el 'damaged record' (bits que ya nos regalan)...")
    tag_high, tag_low = pedir_record_dañado(io)
    print(f"    tag_high=0x{tag_high:05x} tag_low=0x{tag_low:05x}")

    print("[*] pidiendo las claves publicas...")
    pubs = pedir_claves_publicas(io)
    print(f"    tenemos {len(pubs)} claves publicas")

    print("[*] pidiendo firmas de las 5 unidades (4 cada una)...")
    firmas_por_unidad = {i: [] for i in range(UNIT_COUNT)}
    contador_firmas = 0
    for unidad in range(UNIT_COUNT):
        for j in range(SIGN_LIMIT):
            pos = contador_firmas  # tiene que coincidir con el contador del server
            chunk_a = predecir_chunk_a(rc, pos)
            msg = f"firma-{contador_firmas}".encode()
            z, r, s = pedir_firma(io, unidad, msg)
            A_i = (chunk_a << (256 - B_SIZE)) % N
            firmas_por_unidad[unidad].append((r, s, z, A_i))
            contador_firmas += 1
            print(f"    unidad {unidad}, firma {j+1}/{SIGN_LIMIT} (total {contador_firmas})")

    print("[*] a ver si sacamos la clave con el LLL...")
    clave_encontrada = None
    unidad_buena = None
    for unidad, firmas in firmas_por_unidad.items():
        pub = pubs.get(unidad)
        if pub is None:
            continue
        for d in buscar_clave_unidad(firmas):
            # comprobacion rapida con los bits que ya sabiamos
            if (d >> (D_SIZE + A_LOW_SIZE)) != tag_high:
                continue
            if (d & ((1 << A_LOW_SIZE) - 1)) != tag_low:
                continue
            if mult_punto(d, (Gx, Gy)) == pub:
                clave_encontrada = d
                unidad_buena = unidad
                break
        if clave_encontrada is not None:
            break

    if clave_encontrada is None:
        print("[-] pues no ha salido... algo se nos escapa")
        print("    firmas recogidas:", firmas_por_unidad)
        io.close()
        sys.exit(1)

    print(f"[+] clave encontrada, unidad {unidad_buena}: {clave_encontrada}")
    print("[*] mandando el codigo...")
    resultado = mandar_codigo(io, clave_encontrada)
    print(resultado)


if __name__ == "__main__":
    main()
