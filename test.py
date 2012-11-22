from __future__ import division
from sympy import *

x, l, z = symbols("x lambda z", real=True, positive=True)
xi = symbols("xi", real=True)
V, W = map(lambda x: Function(x, real=True, positive=True), "VW")

a_2_func = Lambda((x, xi, l, z), xi**2 + l**2 * V(x) + z**2)
a_0_func = Lambda((x, xi, l, z), W(x))
zero_func = Lambda((x, xi, l, z), 0)

@cacheit
def r(n, a_2=Symbol("a_2"), a_0=Symbol("a_0")):
    if n < 2:
        return zero_func
    if n == 2:
        return a_2(x,xi,l,z) ** (-1)

    rec = lambda k: r(n-k, a_2, a_0)
    # Verified (now for real ...)
    return a_2(x, xi, l, z) ** (-1) * (
                  2 * xi * I * rec(1).diff(x)
                + rec(2).diff(x, 2)
                - (a_0(x, xi, l, z) * rec(2) if n > 3 else 0)
                )

def rr(n):
    return r(n, a_2=a_2_func, a_0=a_0_func)

def trace_exp(n):
    expr = expand(rr(n), deep=False)

    summands = None
    if isinstance(expr, Add):
        summands = expr.args
    else:
        summands = [expr]

    a = (a_2_func(x, xi, l, z) - xi**2)**Rational(1,2)
    b = 1
    c = Wild("c", exclude=[xi])
    n = Wild("n")
    m = Wild("m")

    pattern_0 = c / (a**2 + xi**2)**m
    pattern = c * xi**n / (a**2 + xi**2)**m

    result = []

    for s in summands:
        if s == zero_func:
            val = 0
        else:
            match = pattern.matches(s)
            if match is None:
                match = pattern_0.matches(s)
                # Will fail if the function doesn't have the right format
                match[n] = 0

            val = 0
            if match[n] % 2 == 0:
                # even function over R -> *2, odd function over R -> 0
                x_ = match[n] + 1
                y_ = 2 * match[m] - x_
                val = match[c] * a ** (-y_) \
                      * beta(Rational(x_, 2), Rational(y_, 2))

        result.append(val)

    return sum(result) / (2 * pi)


