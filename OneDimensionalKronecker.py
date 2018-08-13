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
    p = list(map(int, input().split()))
    kronecker(p)

    print("There is the decomposition of polynomial P(x):")
    print(ans)




