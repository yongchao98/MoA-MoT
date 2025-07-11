import numpy as np

def analyze_stbc_diversity():
    """
    Analyzes the diversity order of three 2x2 space-time block codes.

    The script calculates the diversity order for each code by examining the determinant
    of the codeword difference matrix. It uses numerical examples to demonstrate
    rank deficiency for codes with diversity order 1.
    """
    # --- Introduction ---
    print("This script analyzes the diversity order of three space-time codes (Sa, Sb, Sc).")
    print("The diversity order equals the minimum rank of the codeword difference matrix, delta_S = S - S'.")
    print("For a 2x2 code, full diversity (order 2) is achieved if delta_S is always full rank (rank 2) for any two distinct codewords.")
    print("If a non-zero delta_S can have rank 1 (i.e., its determinant is 0), the diversity order is 1.\n")
    print("-" * 65)

    # --- Analysis of Code Sa ---
    print("Analysis for Code Sa = [[x1, x2], [x2, x1]]")
    print("The difference matrix is delta_Sa = [[dx1, dx2], [dx2, dx1]].")
    print("The determinant is det(delta_Sa) = dx1^2 - dx2^2.")
    print("This determinant is zero if dx1 = dx2. We can demonstrate this with an example.")
    dx1_a = 2 + 0j
    dx2_a = 2 + 0j
    delta_Sa = np.array([[dx1_a, dx2_a], [dx2_a, dx1_a]])
    rank_a = np.linalg.matrix_rank(delta_Sa)
    print(f"\nLet's test with dx1 = {dx1_a} and dx2 = {dx2_a}:")
    print(f"delta_Sa =\n{delta_Sa}")
    # Show the numbers in the final equation
    print(f"det(delta_Sa) = ({dx1_a})^2 - ({dx2_a})^2 = {dx1_a**2} - {dx2_a**2} = {dx1_a**2 - dx2_a**2}")
    print(f"The rank of this non-zero matrix is {rank_a}.")
    print("\nConclusion: Since a non-zero difference matrix can have rank 1, the diversity order for Sa is 1.")
    diversity_a = 1
    print("-" * 65)

    # --- Analysis of Code Sb ---
    print("Analysis for Code Sb = [[x1, x2], [x2, x1*]]")
    print("The difference matrix is delta_Sb = [[dx1, dx2], [dx2, dx1*]].")
    print("The determinant is det(delta_Sb) = |dx1|^2 - dx2^2.")
    print("This can be zero if dx2 is a real number and |dx1|^2 = dx2^2 (e.g., from Pythagorean triples).")
    # Differences in QAM coordinates are even integers. A (6,8,10) triple works.
    dx1_b = 6 + 8j
    dx2_b = 10 + 0j
    delta_Sb = np.array([[dx1_b, dx2_b], [dx2_b, np.conj(dx1_b)]])
    rank_b = np.linalg.matrix_rank(delta_Sb)
    print(f"\nLet's test with dx1 = {dx1_b} and dx2 = {dx2_b}:")
    print(f"delta_Sb =\n{delta_Sb}")
    # Show the numbers in the final equation
    print(f"det(delta_Sb) = |{dx1_b}|^2 - ({dx2_b})^2 = {np.abs(dx1_b)**2} - {dx2_b**2} = {np.abs(dx1_b)**2 - dx2_b**2}")
    print(f"The rank of this non-zero matrix is {rank_b}.")
    print("\nConclusion: Since a non-zero difference matrix can have rank 1, the diversity order for Sb is 1.")
    diversity_b = 1
    print("-" * 65)

    # --- Analysis of Code Sc ---
    print("Analysis for Code Sc = [[-x1*, x2], [-x2*, -x1]]")
    print("The difference matrix is delta_Sc = [[-dx1*, dx2], [-dx2*, -dx1]].")
    print("The determinant is det(delta_Sc) = (-dx1*)*(-dx1) - (dx2)*(-dx2*) = |dx1|^2 + |dx2|^2.")
    print("This is a sum of non-negative terms and is zero ONLY if dx1=0 and dx2=0.")
    print("Therefore, any non-zero difference matrix is non-singular and has full rank.")
    dx1_c = 2 + 2j
    dx2_c = 4 - 2j
    delta_Sc = np.array([[-np.conj(dx1_c), dx2_c], [-np.conj(dx2_c), -dx1_c]])
    rank_c = np.linalg.matrix_rank(delta_Sc)
    det_c_val = np.linalg.det(delta_Sc)
    print(f"\nLet's test with dx1 = {dx1_c} and dx2 = {dx2_c}:")
    print(f"delta_Sc =\n{delta_Sc}")
    # Show the numbers in the final equation
    print(f"det(delta_Sc) = |{dx1_c}|^2 + |{dx2_c}|^2 = {np.abs(dx1_c)**2:.1f} + {np.abs(dx2_c)**2:.1f} = {np.abs(dx1_c)**2 + np.abs(dx2_c)**2:.1f}")
    print(f"The calculated determinant is {det_c_val:.1f}, and the rank is {rank_c}.")
    print("\nConclusion: The rank is always 2 for any non-zero inputs. The diversity order for Sc is 2.")
    diversity_c = 2
    print("-" * 65)

    # --- Final Summary ---
    print("\nSUMMARY OF RESULTS")
    print("\n(a) Diversity Orders:")
    print(f"  - Diversity order of Code Sa: {diversity_a}")
    print(f"  - Diversity order of Code Sb: {diversity_b}")
    print(f"  - Diversity order of Code Sc: {diversity_c}")

    print("\n(b) Maximum 'Directivity' (Diversity):")
    print("Code Sc provides the maximum diversity order, which leads to the most reliable communication link against fading.")
    print(f"The maximum diversity order achieved is {diversity_c}.")
    
    final_answer = f"(a) Diversity order for Sa is {diversity_a}, for Sb is {diversity_b}, and for Sc is {diversity_c}. (b) Code Sc provides the maximum diversity."
    print("\n<<<" + final_answer + ">>>")

if __name__ == '__main__':
    analyze_stbc_diversity()