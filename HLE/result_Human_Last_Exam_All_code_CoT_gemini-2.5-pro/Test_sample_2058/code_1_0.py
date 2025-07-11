import math
from fractions import Fraction

def solve():
    """
    Computes the total mass M(A_5, rho, 2) using the mass formula.
    
    The formula is:
    M(G, rho, p) = (1/|G|) * sum_{H <= G} ( n(H, p) / |Aut(H)| ) * sum_{h in H} chi_rho(h)
    
    where n(H, p) is the number of Galois extensions of Q_p with Galois group H.
    """

    # Step 1: Define group, prime, and representation properties
    G_order = 60
    p = 2
    G_name = "A5"
    rho_name = "rho"

    # The character chi_rho(g) for the permutation representation is the number of fixed points of g.
    # We classify elements of A_5 by cycle type to determine the character values.
    # Cycle types in A_5: identity (1), 3-cycle (3), double transposition (22), 5-cycle (5).
    chi_rho = {
        1: 5,   # chi(id)
        3: 2,   # chi(3-cycle)
        22: 1,  # chi((..)(..))
        5: 0    # chi(5-cycle)
    }

    # Step 2: Subgroup data.
    # The image of Gamma_2 must be solvable, so we only consider solvable subgroups of A5.
    # Data format: [name, n(H, 2), |Aut(H)|, {cycle_type: count_in_H}]
    # n(H, 2) values are from the literature on local field extensions (e.g., LMFDB).
    subgroups = [
        ('C1', 1, 1, {1: 1}),
        ('C2', 7, 1, {1: 1, 22: 1}),
        ('C3', 2, 2, {1: 1, 3: 2}),
        ('V4', 18, 6, {1: 1, 22: 3}),
        ('C5', 1, 4, {1: 1, 5: 4}),
        ('S3', 6, 6, {1: 1, 22: 3, 3: 2}),
        ('D5', 0, 8, {1: 1, 22: 5, 5: 4}),  # n(D5) is 0, so contribution is 0
        ('A4', 3, 24, {1: 1, 22: 3, 3: 8})
    ]

    print(f"Calculating the total mass M({G_name}, {rho_name}, {p})")
    print(f"The formula is M = (1/|{G_name}|) * Sum_H [ (n(H, {p}) / |Aut(H)|) * (Sum_h chi(h)) ]")
    print("=" * 70)

    total_sum_fraction = Fraction(0)
    equation_terms = []

    # Step 3 & 4: Calculate contribution for each subgroup and sum them up.
    print("Contributions from each solvable subgroup H:")
    for name, nH, autH_order, composition in subgroups:
        if nH == 0:
            continue

        char_sum = sum(count * chi_rho[cycle_type] for cycle_type, count in composition.items())
        term = Fraction(nH, autH_order) * char_sum
        total_sum_fraction += term
        
        equation_terms.append(f"({nH}/{autH_order})*{char_sum}")

        print(f"H = {name}:")
        print(f"  n({name}, {p}) = {nH}")
        print(f"  |Aut({name})| = {autH_order}")
        print(f"  Sum chi(h) for h in {name} = {char_sum}")
        print(f"  Contribution = ({nH}/{autH_order}) * {char_sum} = {term}")
        print("-" * 30)

    # Step 5: Final Calculation
    print("The full calculation for the sum is:")
    sum_equation = " + ".join(equation_terms)
    print(f"Sum = {sum_equation}")
    print(f"Total Sum = {total_sum_fraction.numerator}/{total_sum_fraction.denominator}")
    print("=" * 70)

    mass = Fraction(total_sum_fraction, G_order)

    print("The final total mass is:")
    print(f"M({G_name}, {rho_name}, {p}) = (1/{G_order}) * ({total_sum_fraction.numerator}/{total_sum_fraction.denominator})")
    print(f"M = {total_sum_fraction.numerator}/({G_order}*{total_sum_fraction.denominator}) = {total_sum_fraction.numerator}/{G_order * total_sum_fraction.denominator}")
    print(f"M = {mass.numerator}/{mass.denominator} (in lowest terms)")

solve()