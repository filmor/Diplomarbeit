
def recurse(t, k_a, k_b, K_a, K_b):
    if len(t) == 0:
        return [[k_a, k_b, K_a, K_b]]

    d = {
        +1: {
            "A_1": [2+k_a, k_b, 1, 1+k_a+k_b],
            "A_2": [1, k_b, 1, k_b],
            "B_1": [2+k_a, k_b, 2+K_a, K_b],
            "B_2": [1, 1+k_b+k_a, 1, 1+K_a+K_b],
            "C":   [1, 1+K_a+K_b, 2+K_a, K_b]
        },
        -1: {
            "A_1": [-1, 1+K_a+K_b, 1, k_b+k_a-1],
            "A_2": [1, K_b, 1, K_b],
            "B_1": [-1, 1+k_a+k_b, K_a, K_b],
            "B_2": [k_a, k_b, 1, K_a+K_b-1],
            "C": [-1, 1+k_a+k_b, K_a, K_b]
        }
    }

    res = []
    for val in d[1 if t[0] >= 0 else -1].itervalues():
        res += recurse(t[1:], *val)
    return res

def calc(*args):
    head = args[0]
    tail = args[1:]

    if head >= 0:
        arr = recurse(tail, 1, 1, 1, 1)
    else:
        arr = recurse(tail, -1, 1, 1, -1)

    return [
        (a[0] + a[1], a[2] + a[3]) for a in arr
    ]

