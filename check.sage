import diamond

BOUND = 15

# coefficient of [Sym^i C] in the identity
def coefficient(k, g, i):
    L = diamond.lefschetz()
    shift = L**(i * (g-k-1))

    # the coefficient is 1 in this case
    if k == g-2 and i == g-1: return shift

    C = gaussian_binomial(2*g-k-i, k + 1-i)(L) \
        - (L^(g-k-1) + L^(g+2*k-3*i)) * gaussian_binomial(2*g-k-i-4, k-i)(L) \
        - (L^(g-k) + L^(g-i) + L^(g+k-2*i) + L^(3*g-3*k-4) + L^(3*g-2*k-4-i) \
            + L^(3*g-k-2*i-4)) * gaussian_binomial(2*g-k-i-4, k-i-1)(L) \
        - (L^(3 * (g-k-1)) + L^(3*(g-k-1) + 1) + L^(3*g-2*k-i-3) \
            + L^(3*g-2*k-i-2)) * gaussian_binomial(2*g-k-i-4, k-i-2)(L) \
        - L^(4*(g-k)-2) * gaussian_binomial(2*g-k-i-4, k-i-3)(L)

    # we check the expression is effective, despite the minus signs
    assert all(c >= 0 for c in C.polynomial.coefficients())

    return C * shift

# testing identity in K_0(Var/k) as equality of Hodge polynomials
for g in range(2, BOUND):
    for k in range(g-1):
        X = diamond.fano_variety_intersection_quadrics_odd(g, k)

        Y = sum(coefficient(k, g, i) * diamond.symmetric_power(i, g) \
                for i in range(k+2))

        print("Testing identity in K_0(Var/k) for g = {}, k = {}".format(g, k))
        assert X == Y
