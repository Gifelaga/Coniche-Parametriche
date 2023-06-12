#################################################
#                                               #
#   Programma sviluppato da Federico Porzio     #
#   Scuola: Arturo Labriola                     #
#   Classe: III H, a.s. 2022/2023               #
#                                               #
#################################################

# Importo la libreria necessaria
import sympy as sp
from sympy import Eq
from sympy import poly
from sympy import And

# Inizializzo i simboli x, y e k
x, y, k = sp.symbols('x y k')

# Prendo in input l'equazione della conica
eq_str = input('\n\nInserire l\'equazione della conica parametrica in forma canonica (es: x**2+y*(k-1)+k*x-2): ')

# Converto la stringa in un'espressione Sympy
eq = sp.S(eq_str)

# Espando eventuali parentesi
new_eq = sp.expand(eq)

peq = poly(new_eq)

# Applico la scombinazione lineare
new_eq = sp.collect(new_eq, k)

# Calcolo le generatrici
genk = new_eq.coeff(k)
gen = new_eq - genk * k

# Stampo le generatrici
print('\n<---------------------------------------------------------------------------------------------------->\n')
print('1. Generatrici: ')
print(gen)
print(genk)
print('\n<---------------------------------------------------------------------------------------------------->\n')

# Calcolo i punti basi
p_base = sp.solve((genk, gen), (x, y))

# Stampo i punti base ottenuti
print('2. Punti base: ')
print(p_base)
print('\n<---------------------------------------------------------------------------------------------------->\n')

# Costruisco la matrice stabilendo i coefficienti nell'equazione, dunque calcolo il determinante
A = eq.coeff(x**2)
B = eq.coeff(x*y)/2
C = eq.coeff(y**2)

p_xy = eq.coeff(x*y)
p = eq - p_xy*x*y
D = p.coeff(x)/2
E = p.coeff(y)/2
F = eq.subs([(x,0), (y,0)])

mat = sp.Matrix([
    [A,B,D],
    [B,C,E],
    [D,E,F]
])
det = mat.det()

# Stampo la matrice ed il determinante
print('3. Matrice: ' + str(mat))
print('\nDeterminante:')
print(det)

# Scovo quando la conica è degenere e stampo il risultato
sol_det=sp.solve(det,k)
print('\nPer k = ' + str(sol_det) + ' la conica è degenere.\n')
CE = sp.Ne(det, 0)

# MAGARI POI CALCOLO IN COSA DEGENERA (ma per il momento non è strettamente necessario)

# Costruisco la matrice MINI e ne calcolo il determinante
mini_mat = sp.Matrix([
    [A,B],
    [B,C]
])
mini_det = mini_mat.det()
print('Matrice dei termini quadratici: ' + str(mini_mat) + '\n')

# Classifico la conica studiando il determinante (><=0) con k come incognita e stampo il risultato:

# Caso parabola
eq = Eq(mini_det,0)
sol = sp.solve((eq,CE),k)
print('*** Per k ∈ ' + str(sol.as_set()) + ' la conica è una parabola')

# Caso iperbole
diseq = mini_det < 0
sol = sp.solve((diseq,CE),k,dict=True)
print('*** Per k ∈ ' + str(sol.as_set()) + " la conica è un'iperbole")

# Sottocaso iperbole equilatera
eq1 = Eq(A, -C)
eq2 = diseq
solutions = sp.solve((eq1,eq2,CE),k)
print('*** Per k ∈ ' + str(solutions.as_set()) + ' la conica è un\'iperbole equilatera')

# Caso ellisse
diseq = mini_det > 0
sol = sp.solve((diseq,CE),k)
print('*** Per k ∈ ' + str(sol.as_set()) + " la conica è un'ellisse")

# Sottocaso circonferenza
eq1 = Eq(A, C)
eq2 = Eq(B, 0)
eq3 = diseq
solutions = sp.solve((eq1,eq2,eq3,CE),k)
print('*** Per k ∈ ' + str(solutions.as_set()) + ' la conica è una circonferenza')

print('\n<---------------------------------------------------------------------------------------------------->\n')
