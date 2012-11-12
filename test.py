from __future__ import division
from sympy import *

x, l, z = symbols("x lambda z", real=True, positive=True)
xi = symbols("xi", real=True)
V, W = map(lambda x: Function(x, real=True, positive=True), "VW")

a_2_func = Lambda((x, xi, l, z), xi**2 + l**2 * V(x) + z**2)
a_0_func = Lambda((x, xi, l, z), W(x))

@cacheit
def D(f, x, n=1):
    return (-I)**n * f.diff(x, n)

@cacheit
def r(n, a_2=Symbol("a_2"), a_0=Symbol("a_0")):
    if n < 2:
        return Lambda((x, xi, l, z), 0)
    if n == 2:
        return a_2(x,xi,l,z) ** (-1)

    rec = lambda k: r(n-k, a_2, a_0)
    # Verified
    return a_2(x, xi, l, z) ** (-1) * (
                  2 * xi * I * D(rec(1), x)
                + 2 * D(rec(2), x, 2)
                - (a_0(x, xi, l, z) * rec(2) if n > 3 else 0)
                )

def rr(n):
  return r(n, a_2=a_2_func, a_0=a_0_func)

