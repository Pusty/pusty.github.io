# https://inria.hal.science/hal-01388108/document


n = 16
IR = IntegerModRing(2**n)
R.<x> = IR[]

t = [-1 for i in range(n+1)]
for i in range(n+1):
    fac = factorial(i)
    for l in range(n+1):
        if (fac % (1<<l)) == 0:
            t[i] = l
       
dn = 0 
for i in range(n+1):       
    if((n-i-t[i]) > 0):
        print(i, n-i-t[i]) 
        dn = i


P = vector(R, 
[2**(n-i-t[i]) * prod([x-j for j in range(2*i)]) for i in range(dn+1)]
+
[prod([x-j for j in range(2*(dn+1))]) ]
)

            

def minimize(inp):
    c = inp
    lg = -1 # anti stuck protection
    while c.degree() > P[-1].degree() and c.degree() != lg:
        lg = c.degree()
        for i in range(0, dn+2):
            exp = 2**(n-i-t[i])
            if i == dn+1:
                exp = 1
            if ZZ(c[c.degree()]) % exp == 0:
                c = c - (R(ZZ(c[c.degree()]) / exp) * x^(c.degree()-P[i].degree())) * P[i]
                break
    return c
            

def newton_invert(i, f):
    if i == 0: return x
    g = newton_invert(i-1, f)
    gd = derivative(g)
    res = g - (gd * (f(g)-x))
    minres = minimize(res)
    print(i, ":", gd, "#", res, "=>", g, " -# ", f(g))
    return minres



f = 123 + 3*x + 42*x^2 + 42*x^3
g = newton_invert(ceil(sqrt(n))+4, f)

# reduction

print(f)
print(g)
print(f(g)(1337))