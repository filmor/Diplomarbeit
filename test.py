from __future__ import division
from sympy import *

x, l, z = symbols("x lambda z", real=True, positive=True)
xi = symbols("xi", real=True)
V, W = map(lambda x: Function(x, real=True, positive=True), "VW")

# a_2, a_0 for the given operator
a_2_func = Lambda((x, xi, l, z), xi**2 + l**2 * V(x) + z**2)
a_0_func = Lambda((x, xi, l, z), W(x))
zero_func = Lambda((x, xi, l, z), 0)

@cacheit
def r(n, a_2=Symbol("a_2"), a_0=Symbol("a_0")):
    if n == 0:
        return zero_func
    if n == 1:
        return zero_func
    if n == 2:
        return a_2(x,xi,l,z) ** (-1)

    rec = lambda k: r(n-k, a_2, a_0)

    return a_2(x, xi, l, z) ** (-1) * (
                  2 * xi * I * rec(1).diff(x)
                + rec(2).diff(x, 2)
                - rec(2) * a_0(x, xi, l, z)
                )

def rr(n):
    return r(n, a_2=a_2_func, a_0=a_0_func)(x,xi,l,z)

def kernel_exp(n):
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

@cacheit
def shifted_harmonic_sum(n, s=Rational(1,2)):
    if n < s:
        return 0
    return shifted_harmonic_sum(n - 1, s) + Rational(1, n - s)

def log_lambda_integration(expr):
    expr = expr.diff(z) / (-2 * z)
    expr = expr.subs(z, 1)
    expr = expand(expr, False)
    if isinstance(expr, Add):
        summands = expr.args
    else:
        summands = [expr]

    k = Wild("k")
    n = Wild("n")
    p = Wild("p", exclude=[l])

    pattern_0 = p / (1 + k * l**2)**((n+5) / 2)
    pattern = pattern_0 * l**n

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

            if match[n] % 2 == 0:
                n_ = match[n]
                k_ = match[k]
                p_ = match[p]
                val = p_ * Rational(1, (n_ + 1) * (n_ + 3)) \
                    * k_ ** (-(n_ + 1)) \
                    * (shifted_harmonic_sum(n_ / 2) - 1 - 2 * ln(2 * k_))
                result.append(val)
            else:
                raise NotImplementedError

    return sum(result)

from sympy.printing.latex import LatexPrinter
class MyLatexPrinter(LatexPrinter):
    _default_settings = LatexPrinter._default_settings

    def __init__(self, **kwargs):
        super(MyLatexPrinter, self).__init__(kwargs)

    def _print_Derivative(self, expr):
        dim = len(expr.variables)

        if len(set(expr.variables)) != 1 or not isinstance(expr.expr, Function):
            return super(MyLatexPrinter, self)._print_Derivative(expr)

        deriv_expr = ""
        if dim == 1:
            deriv_expr = "'"
        elif dim == 2:
            deriv_expr = "''"
        else:
            deriv_expr = r"^{(%s)}" % dim

        return self._print_Function(expr.expr, deriv=deriv_expr)

    def _print_Function(self, expr, exp=None, deriv=""):
        name = expr.func.__name__
        if deriv is "" and name not in ["W", "V", "a_2", "a_0"]:
            return super(MyLatexPrinter, self)._print_Function(expr)
        else:
            args = ",".join(str(self._print(arg)) for arg in expr.args)
            exp = "^{%s}" % exp if exp else ""
            return "%s%s(%s)%s" % (name, deriv, args, exp)

def format_expr(expr):
    return MyLatexPrinter(mode="plain", fold_func_brackets=True,
            fold_frac_powers=True).doprint(expr)

