import sympy

def solve_stbc_diversity():
    """
    This script analyzes the diversity order of three 2x2 space-time block codes (STBCs)
    and determines which one offers the maximum diversity gain (directivity).

    The diversity order is found by examining the determinant of the codeword difference matrix.
    If the determinant is zero only for zero error vectors, the code achieves full diversity (order 2).
    Otherwise, the diversity order is 1.
    """

    # Define symbolic complex variables for the error symbols e1 and e2
    # e_i represents the difference between two distinct symbols, i.e., e_i = x_i - x'_i
    e1, e2 = sympy.symbols('e1 e2', complex=True)

    print("--- Analysis of Space-Time Block Codes ---")

    # --- Part (a): Diversity Order Analysis ---
    print("\n--- (a) Diversity Order for each code ---")

    # --- Code S_a ---
    print("\n1. Analysis of Code S_a = [[x1, x2], [x2, x1]]")
    # The difference matrix Delta_S_a is [[e1, e2], [e2, e1]]
    det_a = e1**2 - e2**2
    print(f"The determinant of the difference matrix is: det(Delta_S_a) = {det_a}")
    print("To find the diversity order, we check if the determinant can be zero for non-zero error symbols.")
    print(f"The equation is: {det_a} = 0, which means e1**2 = e2**2.")
    print("This condition is met if e1 = e2 or e1 = -e2. We can easily find non-zero error symbols that satisfy this.")
    print("For example, consider the non-zero error vector e = [e1, e2] = [1, 1].")
    print(f"Plugging these values into the determinant equation: (1)**2 - (1)**2 = 1 - 1 = 0.")
    print("Since the determinant can be zero for non-zero errors, the code does not have full diversity.")
    print("Result: The diversity order for S_a is 1.")

    # --- Code S_b ---
    print("\n2. Analysis of Code S_b = [[x1, x2], [x2, x1*]]")
    # The difference matrix Delta_S_b is [[e1, e2], [e2, conjugate(e1)]]
    # The determinant is e1*conjugate(e1) - e2*e2
    det_b_expr_str = "|e1|**2 - e2**2"
    print(f"The determinant of the difference matrix is: det(Delta_S_b) = {det_b_expr_str}")
    print("We check if the determinant can be zero for non-zero error symbols.")
    print(f"The equation is: {det_b_expr_str} = 0, which means |e1|**2 = e2**2.")
    print("For this to hold, e2**2 must be a non-negative real number. This implies e2 must be a real number.")
    print("Let's see if we can find a non-zero complex e1 and a non-zero real e2 that satisfy this.")
    print("Consider the non-zero error vector where e1 = 3 + 4j and e2 = 5.")
    e1_val = 3 + 4j
    e2_val = 5
    det_b_val = abs(e1_val)**2 - e2_val**2
    print(f"Plugging these values into the determinant equation: |{e1_val}|**2 - ({e2_val})**2 = {abs(e1_val)**2} - {e2_val**2} = {det_b_val}.")
    print("Since the determinant can be zero for non-zero errors, the code does not have full diversity.")
    print("Result: The diversity order for S_b is 1.")

    # --- Code S_c ---
    print("\n3. Analysis of Code S_c = [[-x1*, x2], [-x2*, -x1]]")
    # The difference matrix Delta_S_c is [[-e1*, e2], [-e2*, -e1]]
    # The determinant is (-e1*)*(-e1) - e2*(-e2*) = |e1|**2 + |e2|**2
    det_c_expr_str = "|e1|**2 + |e2|**2"
    print(f"The determinant of the difference matrix is: det(Delta_S_c) = {det_c_expr_str}")
    print("We check if the determinant can be zero for non-zero error symbols.")
    print(f"The equation is: {det_c_expr_str} = 0.")
    print("In this final equation, the numbers are coefficients 1 for |e1|**2, 1 for |e2|**2, and 0 on the right side.")
    print("Since |e1|**2 and |e2|**2 are both non-negative, their sum can only be zero if both terms are zero.")
    print("  |e1|**2 = 0 implies e1 = 0.")
    print("  |e2|**2 = 0 implies e2 = 0.")
    print("The determinant is zero only for the zero error vector [e1, e2] = [0, 0].")
    print("For any non-zero error vector, the determinant is strictly positive, so the code has full diversity.")
    print("Result: The diversity order for S_c is 2.")

    # --- Part (b): Maximum Directivity ---
    print("\n--- (b) Code with Maximum Directivity ---")
    print("\n'Directivity' is interpreted as diversity gain. Maximum directivity corresponds to the highest diversity order.")
    print("The maximum possible diversity order for a 2x1 MISO system is 2.")
    print("\nSummary of diversity orders:")
    print("  - Diversity order of S_a: 1")
    print("  - Diversity order of S_b: 1")
    print("  - Diversity order of S_c: 2")
    print("\nConclusion: Code S_c achieves the maximum diversity order of 2 and therefore provides the maximum directivity.")

if __name__ == '__main__':
    solve_stbc_diversity()