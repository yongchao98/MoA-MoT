import math

def combinations(n, k):
    """Computes the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    # Use math.comb for Python 3.8+ for efficiency and accuracy
    if hasattr(math, 'comb'):
        return math.comb(n, k)
    # Fallback for older versions
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve_lattice_problems():
    """
    Solves the three parts of the lattice theory problem and prints the results.
    """
    # Part (a) Calculation
    n_a = 16
    c16_4 = combinations(n_a, 4)
    c16_8 = combinations(n_a, 8)
    c16_12 = combinations(n_a, 12)
    c16_16 = combinations(n_a, 16)
    ans_a = c16_4 + c16_8 + c16_12 + c16_16
    
    # Part (b) is a logical deduction
    ans_b = "no"

    # Part (c) Calculation
    n_c = 12
    c12_4 = combinations(n_c, 4)
    c12_8 = combinations(n_c, 8)
    c12_12 = combinations(n_c, 12)
    num_kernels = c12_4 + c12_8 + c12_12
    ans_c = 2 * num_kernels

    # The final output needs to follow the requested format.
    # The intermediate steps are shown for clarity as requested.
    print("Step-by-step calculations:")
    print("\nPart (a): Number of 2-isotropic vectors for Z^16")
    print(f"The number of vectors is the sum of binomial coefficients C(16, k) for k in {{4, 8, 12, 16}}.")
    print(f"C(16, 4) = {c16_4}")
    print(f"C(16, 8) = {c16_8}")
    print(f"C(16, 12) = C(16, 4) = {c16_12}")
    print(f"C(16, 16) = {c16_16}")
    print(f"The final sum is: {c16_4} + {c16_8} + {c16_12} + {c16_16} = {ans_a}")

    print("\nPart (b): Possibility of an even 3-neighbor for Z^8")
    print("The reasoning is based on a proof by contradiction using lattice indices.")
    print(f"The answer is '{ans_b}'.")
    
    print("\nPart (c): Number of unimodular 2-neighbors of Z^12")
    print(f"The number of defining kernels is C(12, k) for k in {{4, 8, 12}}.")
    print(f"C(12, 4) = {c12_4}")
    print(f"C(12, 8) = C(12, 4) = {c12_8}")
    print(f"C(12, 12) = {c12_12}")
    print(f"Number of kernels = {c12_4} + {c12_8} + {c12_12} = {num_kernels}")
    print(f"Each kernel corresponds to 2 neighbors, so the total number is 2 * {num_kernels} = {ans_c}")

    print("\n---")
    print("Final Answer:")
    # Printing the final answer in the format (a) [Numerical]; (b) [yes/no]; (c) [numerical].
    final_answer_string = f"(a) [{ans_a}]; (b) [{ans_b}]; (c) [{ans_c}]."
    print(final_answer_string)
    
if __name__ == '__main__':
    solve_lattice_problems()
