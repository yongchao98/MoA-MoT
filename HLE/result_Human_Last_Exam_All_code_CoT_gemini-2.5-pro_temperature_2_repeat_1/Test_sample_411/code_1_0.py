import numpy as np

def solve_stbc_diversity():
    """
    Analyzes the diversity order of three STBCs and identifies the one with maximum directivity.
    """
    print("This script analyzes the diversity order for three space-time block codes S_a, S_b, and S_c.")

    # Part (a): Diversity Order Analysis
    print("\n" + "="*60)
    print("Part (a): Analysis of Diversity Order")
    print("="*60)
    print("The diversity order of a space-time block code (STBC) is determined by the minimum")
    print("rank of the codeword difference matrix, Delta_S, over all pairs of distinct codewords.")
    print("For the given 2x2 STBCs, full diversity (order 2) is achieved if and only if")
    print("the determinant of Delta_S is non-zero for any non-zero error vector e = [e1, e2].")
    print("\nWe analyze this for each code and use a QPSK constellation as a concrete example.")

    # Define QPSK constellation and the set of difference symbols (error symbols)
    C = [1+1j, 1-1j, -1+1j, -1-1j]
    E_set = set()
    for c1 in C:
        for c2 in C:
            E_set.add(c1 - c2)
    # The error alphabet includes 0, which we need for e1 or e2 being zero.
    E_full = sorted(list(E_set), key=lambda x: (abs(x), np.angle(x))) 

    print("\nQPSK constellation points:", C)
    print("Possible error symbols (e_k):", E_full)

    # --- Analysis for Code S_a ---
    print("\n--- Analysis for Code S_a ---")
    print("S_a = [ [x1, x2], [x2, x1] ]")
    print("The difference matrix is Delta_S_a = [ [e1, e2], [e2, e1] ]")
    print("The determinant is det(Delta_S_a) = e1**2 - e2**2.")
    print("This determinant is zero if e1 = e2 or e1 = -e2. We can easily find non-zero error symbols that satisfy this.")
    
    e1_a = 2+0j
    e2_a = 2+0j
    det_a = e1_a**2 - e2_a**2
    delta_S_a = np.array([[e1_a, e2_a], [e2_a, e1_a]])
    rank_a = np.linalg.matrix_rank(delta_S_a)

    print("\nExample demonstrating rank deficiency for S_a:")
    print(f"  Let e1 = {e1_a} and e2 = {e2_a}. This is a non-zero error vector.")
    print(f"  The determinant equation is: ({e1_a})^2 - ({e2_a})^2 = {e1_a**2} - {e2_a**2} = {det_a}")
    print(f"  Since the determinant is 0, the matrix is singular and its rank is {rank_a}.")
    print("\nConclusion: The minimum rank is 1. Diversity Order for S_a is 1.")

    # --- Analysis for Code S_b ---
    print("\n--- Analysis for Code S_b ---")
    print("S_b = [ [x1, x2], [x2, x1*] ]")
    print("The difference matrix is Delta_S_b = [ [e1, e2], [e2, e1*] ], where e1* is the conjugate of e1.")
    print("The determinant is det(Delta_S_b) = e1*e1* - e2**2 = |e1|**2 - e2**2.")
    print("This can be zero if we find a non-zero pair (e1, e2) where the squared magnitude of e1 equals the square of e2.")
    
    e1_b = 2j
    e2_b = 2+0j
    det_b = abs(e1_b)**2 - e2_b**2
    delta_S_b = np.array([[e1_b, e2_b], [e2_b, np.conj(e1_b)]])
    rank_b = np.linalg.matrix_rank(delta_S_b)
    
    print("\nExample demonstrating rank deficiency for S_b:")
    print(f"  Let e1 = {e1_b} and e2 = {e2_b}. Both are valid non-zero error symbols from the QPSK difference set.")
    print(f"  The determinant equation is: |{e1_b}|^2 - ({e2_b})^2 = {abs(e1_b)**2:.1f} - {e2_b**2} = {det_b.real:.1f}")
    print(f"  Since the determinant is 0, the matrix is singular and its rank is {rank_b}.")
    print("\nConclusion: The minimum rank is 1. Diversity Order for S_b is 1.")
    
    # --- Analysis for Code S_c ---
    print("\n--- Analysis for Code S_c ---")
    print("S_c = [ [-x1*, x2], [-x2*, -x1] ]")
    print("The difference matrix is Delta_S_c = [ [-e1*, e2], [-e2*, -e1] ].")
    print("The determinant is det(Delta_S_c) = (-e1*)*(-e1) - (e2)*(-e2*) = |e1|**2 + |e2|**2.")
    print("This expression is the sum of two non-negative numbers. It is zero if and only if |e1|^2 = 0 and |e2|^2 = 0, which means e1 = 0 and e2 = 0.")
    print("For any non-zero error vector e = [e1, e2], the determinant is always strictly positive.")
    
    e1_c = 2+0j
    e2_c = 2j
    det_c = abs(e1_c)**2 + abs(e2_c)**2
    delta_S_c = np.array([[-np.conj(e1_c), e2_c], [-np.conj(e2_c), -e1_c]])
    rank_c = np.linalg.matrix_rank(delta_S_c)

    print("\nExample calculation for S_c with a non-zero error vector:")
    print(f"  Let e1 = {e1_c} and e2 = {e2_c}.")
    print(f"  The determinant equation is: |{e1_c}|^2 + |{e2_c}|^2 = {abs(e1_c)**2:.1f} + {abs(e2_c)**2:.1f} = {det_c:.1f}")
    print(f"  The determinant is non-zero, so the matrix is full rank ({rank_c}).")
    print("\nConclusion: The minimum rank is 2. Diversity Order for S_c is 2 (Full Diversity).")

    # Final summary for part (a)
    print("\n" + "*"*60)
    print("Summary for Part (a): Diversity Orders")
    print("  - Diversity order for S_a: 1")
    print("  - Diversity order for S_b: 1")
    print("  - Diversity order for S_c: 2")
    print("*"*60)
    
    # Part (b): Maximum Directivity
    print("\n" + "="*60)
    print("Part (b): Which code provides the maximum directivity?")
    print("="*60)
    print("In the context of STBCs, 'directivity' or antenna gain is achieved through diversity.")
    print("A higher diversity order results in better performance against channel fading.")
    print("Comparing the diversity orders, code S_c (with diversity 2) is superior to S_a and S_b (both with diversity 1).")
    print("\nTherefore, code S_c provides the maximum directivity.")

solve_stbc_diversity()