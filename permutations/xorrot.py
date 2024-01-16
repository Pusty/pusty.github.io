# https://marc-b-reynolds.github.io/math/2017/10/13/XorRotate.html
# The invertibility of the XOR of rotations of a binary word - https://doi.org/10.1080/00207161003596708

from sympy.matrices import Matrix, eye, zeros, ones, diag, GramSchmidt
from sympy import *
from sympy.solvers.solveset import linsolve

BITS = 4

OFFSET = (1, 3) # I + C^a + C^b

# can also replace this with rol
def ror(n,rotations,width=BITS):
    return (2**width-1)&(n>>rotations|n<<(width-rotations))


def applyMatrix(matrix, value):
    out = 0
    for y in range(BITS):
        o = 0
        for x in range(BITS):
            b = (value>>x)&1
            o = (o+(matrix[y, x] * b))%2
        out = out | (int(o)<<y)
    return out
    
def applyDirect(value):
    return value ^ ror(value, OFFSET[0]) ^ ror(value, OFFSET[1])

def bits(value):
    out = []
    for x in range(BITS):
        out.append((value>>x)&1)
    return Matrix(out)
    
def debit(matrix):
    out = 0
    for x in range(BITS):
        out = out | (int(matrix[x])<<x)
    return out
    

def applyMatrix2(matrix, value):
    return debit((matrix@bits(value))%2)
    

def invert(M):
    return M.inv_mod(2)

col = []
for y in range(BITS):
    row = [0 for _ in range(BITS)]
    for x in range(BITS):
        if y == x:
            row[x] = row[x] ^ 1
        if (y+OFFSET[0])%BITS == x:
            row[x] = row[x] ^ 1
        if (y+OFFSET[1])%BITS == x:
            row[x] = row[x] ^ 1
    col.append(row)
    
M = Matrix(col)

pprint(M)

# if nullspace is trivial linear transform is injective
# since we map Z^n_2 -> Z^n_2 this is enough for proof of being bijective -> is permutation
# print(M.nullspace())

print(applyDirect(13+42), applyDirect(13) + applyDirect(42))
print(applyDirect(applyDirect(42)))

print("=====")

Minv = invert(M)

pprint(Minv)
print("=====")
