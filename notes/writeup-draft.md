# WRITEUP

FAILS FROM THE CHALLENGE:
1. RNG inseguro
2. Fallo de datos: Al darle a la opción 5 del menu (data panel) 9 veces se reconstruye el estado de la librería de randcrack ya que se necesitan solo 624 salidas y al hacerlo 9 veces nos da 702 salidas.
   
ATAQUE:
1. Miré los puntos que podrían dar error en el código (Los de arriba)
2. Había que hacer un script con una serie de pasos para que interactuara con la ip, asi que se abre el puerto y en otro archvio hice el script. 
3. El script fue hacer que hiciese 9 nueve veces la opcion 5
4. Predecir un trozo del nonce
5. Cada firma tiene r,s,z (que te dan) y el nonce, se aplica formula para encontrar t y c
6. Se hace una matriz con la libreria de LLL de python para reducir de tamaño las operaciones
7. Se enuentra el otro trozo de nonce
8. Se verifica con la clave pública 
9.  Se ejecuta el script abriendo la ip en consola y enseñando el panel del código que nos dan, se selecciona la opción 6 de mandar respuesta y se envía y aparece la flag 
    
## Resultado

Ejecutando el solver contra el reto real:

```
$ sage -python solve.py
[*] Collecting data-panel leaks (9x)...
    read batch 1/9, total entries: 78
    read batch 2/9, total entries: 156
    read batch 3/9, total entries: 234
    read batch 4/9, total entries: 312
    read batch 5/9, total entries: 390
    read batch 6/9, total entries: 468
    read batch 7/9, total entries: 546
    read batch 8/9, total entries: 624
    read batch 9/9, total entries: 702
[*] Cracking MT19937 state...
    alignment check: OK
[*] Reading damaged record (sanity-check bits of every secret)...
    tag_high=0x87fdf tag_low=0x081f7
[*] Reading public keys...
    got 5 public keys
[*] Requesting signatures for all units...
    unit 0 sig 1/4 (total 1)
    unit 0 sig 2/4 (total 2)
    unit 0 sig 3/4 (total 3)
    unit 0 sig 4/4 (total 4)
    unit 1 sig 1/4 (total 5)
    unit 1 sig 2/4 (total 6)
    unit 1 sig 3/4 (total 7)
    unit 1 sig 4/4 (total 8)
    unit 2 sig 1/4 (total 9)
    unit 2 sig 2/4 (total 10)
    unit 2 sig 3/4 (total 11)
    unit 2 sig 4/4 (total 12)
    unit 3 sig 1/4 (total 13)
    unit 3 sig 2/4 (total 14)
    unit 3 sig 3/4 (total 15)
    unit 3 sig 4/4 (total 16)
    unit 4 sig 1/4 (total 17)
    unit 4 sig 2/4 (total 18)
    unit 4 sig 3/4 (total 19)
    unit 4 sig 4/4 (total 20)
[*] Running HNP lattice attack per unit...
[+] Recovered secret for unit 0: 61510954339560518939547126329213547641079480598258146797275258950427988754935
[*] Submitting...

  Accepted for unit #0.
```
```
COMPFEST18{b1as3d_n0nc3_mt_r3c0v3ry_lll_hnp_go_brr_727e3a9724b244c1}
```
