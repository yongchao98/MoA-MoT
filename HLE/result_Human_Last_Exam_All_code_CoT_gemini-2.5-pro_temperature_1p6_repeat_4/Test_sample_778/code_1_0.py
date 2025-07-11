import fractions

def solve_problem():
    """
    Calculates the groupoid cardinality of the inserter Ins(F,G).

    The problem asks for the groupoid cardinality of the inserter of two functors
    F, G: BZ/nZ -> BZ/mZ, where the functors are induced by the maps
    f(x) = ax and g(x) = bx.

    The objects of the inserter category Ins(F,G) can be identified with the
    elements of Z/mZ. A morphism exists from an object k1 to an object k2
    if there is a u in Z/nZ such that k2 - k1 = (b-a)u mod m.

    The structure of this groupoid is a disjoint union of connected components.
    The groupoid cardinality is given by the sum over these components of the
    reciprocal of the size of the automorphism group of any object in that
    component.

    A detailed derivation shows that this value simplifies to the fraction m/n,
    provided that F and G are well-defined functors (i.e., an = 0 mod m and
    bn = 0 mod m), which is true for the given numbers.
    The specific values of a and b do not affect the final result.
    """
    n = 54669191328000
    m = 1470432000
    a = 991
    b = 223
    
    # The groupoid cardinality is given by the fraction m / n.
    cardinality_fraction = fractions.Fraction(m, n)
    
    # The problem asks to output each number in the final equation.
    # We will print the original fraction and its simplified form.
    print(f"{m} / {n} = {cardinality_fraction.numerator} / {cardinality_fraction.denominator}")

solve_problem()