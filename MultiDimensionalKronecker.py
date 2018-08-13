def gorner(p, x):
    sum = p[0]
    for i in range(1, len(p)):
        sum = x * sum + p[i]
    return sum


def find_dividers(x):
    dividers = []
    if x == 0:
        dividers.append(0)
    else:
        for y in range(1, abs(x) // 2 + 1):
            if x % y == 0:
                dividers.append(y)
                dividers.append(-y)
        dividers.append(x)
        dividers.append(-x)

    return dividers


def check(p):
    flag = True
    for x in p:
        from math import ceil
        if x != ceil(x):
            flag = False
            break
    return flag


ans = []


def kronecker(p):
    n = len(p) - 1
    m = n // 2 + 1
    values = []

    for x in range(0, m):
        values.append(x)
        if len(values) < m and x != 0:
            values.append(-x)
        if len(values) == m:
            break

    intrpl = {}
    max_size = 1
    for x in values:
        y = gorner(p, x)
        dividers = find_dividers(y)
        intrpl[x] = dividers
        max_size *= len(dividers)

    mas_a = []
    for x in intrpl.keys():
        mas_a.append(x)
    uniq_a_all = []
    start = 3
    if m == 2:
        start = 2
    for L in range(start, len(mas_a) + 1):
        import itertools
        for subset in itertools.combinations(mas_a, L):
            uniq_a_all.append(list(subset))
    flag_all = False
    for a in uniq_a_all:
        uniq_b = []
        while True:
            b = []
            for x in a:
                length = len(intrpl[x])
                import random
                idx = random.randint(0, length - 1)
                b.append(intrpl[x][idx])
            if b in uniq_b:
                continue
            else:
                uniq_b.append(b)
                from scipy import interpolate
                q_inter = interpolate.lagrange(a, b)
                q = []
                for x in q_inter.coefficients:
                    q.append(int(x))
                from numpy import polydiv
                quot, rem = polydiv(p, q)
                if len(rem) == 1 and rem == [0.]:
                    if len(q) == len(p) or len(quot) == len(p):
                        continue
                    flag = check(q)
                    if not flag:
                        continue
                    p1 = []
                    for x in quot:
                        p1.append(float(x))
                    flag = check(p1)
                    if not flag:
                        continue
                    for i in range(len(p1)):
                        p1[i] = int(p1[i])
                    if len(q) > 2:
                        kronecker(q)
                    else:
                        ans.append(q)
                    if len(p1) > 2:
                        kronecker(p1)
                    else:
                        ans.append(p1)
                    flag_all = True
                    break
            if len(uniq_b) == max_size:
                break
        if flag_all:
            break
    if not flag_all:
        ans.append(p)


if __name__ == '__main__':
    print('Enter the polynomial: ', end='')

    from sympy import Poly, Symbol

    polynom = Poly(input())

    d = polynom.degree() + 1
    a = Symbol('a')

    symbols = []
    for symbol in polynom.free_symbols:
        symbols.append(str(symbol))

    symbols.sort()

    replace_symbols = {}
    j = 0
    for symbol in symbols:
        replace_symbols[symbol] = a ** (d ** j)
        j += 1

    for symbol in replace_symbols:
        polynom = polynom.subs(symbol, replace_symbols[symbol])

    p_new = (polynom._sorted_args)[0]._sorted_args
    f = 0
    for i in range(len(p_new) - 1, -1, -1):
        if i == 0 and str(p_new[i]).find('a') == -1:
            f += int(p_new[i])
        else:
            f += Poly(str(p_new[i]))
    p_coef = f.all_coeffs()
    for i in range(len(p_coef)):
        p_coef[i] = int(p_coef[i])

    kronecker(p_coef)

    ans.sort(key=len)
    for i in range(len(ans)):
        deg = len(ans[i]) - 1
        x = 1
        while deg > x:
            x *= d
        flag = False
        if 0.5 < deg / x < 1:
            flag = True
        if flag:
            for j in range(0, i):
                if len(ans[j]) - 1 + deg == x:
                    from numpy import polymul
                    mult = polymul(ans[i], ans[j])
                    new_coef = []
                    for x in mult:
                        new_coef.append(int(x))
                    ans[i] = new_coef
                    ans.remove(ans[j])
                    break

    final = []
    for polynomial in ans:
        q = Poly(polynomial, a)
        for symbol in reversed(symbols):
            q = q.subs(replace_symbols[symbol], symbol)
        final.append(q)

    print(final)
