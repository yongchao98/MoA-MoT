import math

def compute_mass():
    """
    Computes the total mass M(A_5, rho, 2).
    """

    # Properties of subgroups of A_4
    # N: Number of Galois extensions of Q_2 with this Galois group. Source: Jones-Roberts Database / LMFDB.
    # n: Number of subgroups of A_4 isomorphic to this group.
    # aut_order: Order of the automorphism group of the subgroup.
    subgroups_data = [
        {'name': 'Trivial {e}', 'N': 1, 'n': 1, 'aut_order': 1},
        {'name': 'C_2', 'N': 7, 'n': 3, 'aut_order': 1},
        {'name': 'C_3', 'N': 1, 'n': 4, 'aut_order': 2},
        {'name': 'V (Klein-4)', 'N': 19, 'n': 1, 'aut_order': 6}, # Aut(V) is S_3
        {'name': 'A_4', 'N': 2, 'n': 1, 'aut_order': 24}  # Aut(A_4) is S_4
    ]

    print("Step 1: The problem reduces to computing M(A_4, triv, 2) = |Hom(Gamma_2, A_4)| / |A_4|.")
    print("Step 2: Calculate |Hom(Gamma_2, A_4)| by summing contributions from each subgroup type.")
    print("-" * 30)

    total_homomorphisms = 0
    for H in subgroups_data:
        name = H['name']
        N = H['N']
        n = H['n']
        aut_order = H['aut_order']
        
        num_injections = n * aut_order
        contribution = N * num_injections
        total_homomorphisms += contribution
        
        print(f"Subgroup type: {name}")
        print(f"  Number of Q_2 extensions N(H) = {N}")
        print(f"  Number of injections |Inj(H, A_4)| = n(H) * |Aut(H)| = {n} * {aut_order} = {num_injections}")
        print(f"  Contribution to |Hom(Gamma_2, A_4)| = N(H) * |Inj(H, A_4)| = {N} * {num_injections} = {contribution}")
        print("-" * 30)

    print(f"Total number of homomorphisms |Hom(Gamma_2, A_4)| = {total_homomorphisms}")

    order_A4 = 12
    print(f"The order of A_4 is {order_A4}.")
    
    mass_A4 = total_homomorphisms / order_A4
    
    print(f"M(A_4, triv, 2) = {total_homomorphisms} / {order_A4} = {int(mass_A4)}")

    mass_A5 = mass_A4
    print("\nStep 3: M(A_5, rho, 2) = M(A_4, triv, 2) due to the property of induced representations.")
    print(f"So, M(A_5, rho, 2) = {int(mass_A5)}")

    # Express as a fraction in lowest terms
    numerator = int(mass_A5)
    denominator = 1
    
    # In this case, gcd is trivial, but it's good practice
    common_divisor = math.gcd(numerator, denominator)
    num_simple = numerator // common_divisor
    den_simple = denominator // common_divisor
    
    print(f"\nThe final answer as a fraction in lowest terms is {num_simple}/{den_simple}.")

compute_mass()
<<<16/1>>>