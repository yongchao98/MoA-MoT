import sympy

def analyze_stbc():
    """
    Analyzes the diversity order of three space-time block codes
    and determines which one offers the best performance.
    """
    # Define complex symbols for the QAM constellation points and their differences (error symbols)
    e1, e2 = sympy.symbols('e_1 e_2', complex=True, nonzero=True)
    
    print("--- Analysis of Space-Time Block Codes ---")
    print("The diversity order is determined by the rank of the codeword difference matrix E.")
    print("A code has full diversity (order 2) if det(E) is never zero for distinct codewords.\n")

    # --- Code Sa ---
    print("--- (1) Analysis for Code Sa ---")
    # E_a = [[e1, e2], [e2, e1]]
    det_Ea = e1**2 - e2**2
    print(f"The determinant of the difference matrix E_a is: det(E_a) = {det_Ea}")
    print("The equation is (1)*e_1**2 + (-1)*e_2**2.")
    print("This determinant can be zero if e_1 = e_2 or e_1 = -e_2.")
    print("Since e_1 and e_2 are differences of M-QAM symbols, it's possible to find distinct codewords where this condition holds.")
    print("For example, using 4-QAM, if x_1=1+j, x_2=1+j and x'_1=-1-j, x'_2=-1-j, then e_1 = 2+2j and e_2 = 2+2j.")
    print("Since det(E_a) can be zero, the code is not full diversity.")
    print("Diversity Order for S_a: 1\n")
    diversity_a = 1

    # --- Code Sb ---
    print("--- (2) Analysis for Code Sb ---")
    # E_b = [[e1, e2], [e2, conj(e1)]]
    det_Eb = sympy.Abs(e1)**2 - e2**2
    print(f"The determinant of the difference matrix E_b is: det(E_b) = |e_1|**2 - e_2**2")
    print("The equation is (1)*|e_1|**2 + (-1)*e_2**2.")
    print("This determinant can be zero if |e_1|**2 = e_2**2.")
    print("This is possible. For instance, if e_2 is a real number (e.g., e_2=2, from x_2=1+j, x'_2=-1+j),")
    print("and e_1 is chosen such that |e_1| = |e_2| (e.g., e_1=2, from x_1=1+j, x'_1=-1+j).")
    print("Since det(E_b) can be zero, the code is not full diversity.")
    print("Diversity Order for S_b: 1\n")
    diversity_b = 1

    # --- Code Sc ---
    print("--- (3) Analysis for Code Sc ---")
    # E_c = [[-conj(e1), e2], [-conj(e2), -e1]]
    det_Ec = sympy.Abs(e1)**2 + sympy.Abs(e2)**2
    print(f"The determinant of the difference matrix E_c is: det(E_c) = |e_1|**2 + |e_2|**2")
    print("The equation is (1)*|e_1|**2 + (1)*|e_2|**2.")
    print("Since e_1 and e_2 cannot both be zero (as codewords are distinct), |e_1|**2 + |e_2|**2 is always greater than zero.")
    print("The determinant is never zero for any pair of distinct codewords.")
    print("Therefore, the code is full-rank and achieves full diversity.")
    print("Diversity Order for S_c: 2\n")
    diversity_c = 2

    # --- Summary and Conclusion ---
    print("--- (a) Summary of Diversity Orders ---")
    print(f"Diversity Order of S_a: {diversity_a}")
    print(f"Diversity Order of S_b: {diversity_b}")
    print(f"Diversity Order of S_c: {diversity_c}\n")

    print("--- (b) Code with Maximum Directivity ---")
    print("Maximum directivity, in the context of space-time codes, refers to the best performance.")
    print("Performance is primarily determined by the diversity order.")
    print(f"Code S_c has the highest diversity order ({diversity_c}), while S_a and S_b have lower diversity orders ({diversity_a} and {diversity_b}).")
    print("Therefore, S_c provides the maximum directivity (best performance).")

if __name__ == '__main__':
    analyze_stbc()