import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve_part_a():
    """
    Calculates the number of 2-isotropic vectors for Z^16.
    This corresponds to counting non-zero vectors in (Z/2Z)^16
    with weight divisible by 4.
    """
    n = 16
    total_vectors = 0
    print("Part (a) calculation:")
    print("Counting non-zero vectors v in (Z/2Z)^16 with weight k divisible by 4.")
    print("Number of vectors = C(16, 4) + C(16, 8) + C(16, 12) + C(16, 16)")
    
    terms = []
    # k must be a non-zero multiple of 4
    for k in range(4, n + 1, 4):
        term = combinations(n, k)
        terms.append(f"C(16, {k}) = {term}")
        total_vectors += term
    
    print(" + ".join(terms))
    print(f"= {' + '.join([str(combinations(n, k)) for k in range(4, n + 1, 4)])} = {total_vectors}")
    return total_vectors

def solve_part_b():
    """
    Determines if a 3-neighbor of Z^8 can be even.
    Based on logical deduction, it is not possible.
    """
    print("\nPart (b) reasoning:")
    print("For a neighbor N_3(x) of Z^8 to be even, x must satisfy conditions that lead to a contradiction:")
    print("1. All components x_i of x must not be divisible by 3.")
    print("2. The sum of squares x.x must be divisible by 3.")
    print("Condition (1) implies x.x = 8 (mod 3), which contradicts condition (2).")
    print("Therefore, it is not possible.")
    return "no"

def solve_part_c():
    """
    Calculates the number of unimodular 2-neighbors of Z^12.
    This corresponds to counting non-zero vectors in (Z/2Z)^12
    with weight divisible by 4.
    """
    n = 12
    total_neighbors = 0
    print("\nPart (c) calculation:")
    print("Counting non-zero vectors v in (Z/2Z)^12 with weight k divisible by 4.")
    print("Number of neighbors = C(12, 4) + C(12, 8) + C(12, 12)")
    
    terms = []
    # k must be a non-zero multiple of 4
    for k in range(4, n + 1, 4):
        term = combinations(n, k)
        terms.append(f"C(12, {k}) = {term}")
        total_neighbors += term

    print(" + ".join(terms))
    print(f"= {' + '.join([str(combinations(n, k)) for k in range(4, n + 1, 4)])} = {total_neighbors}")
    return total_neighbors

# --- Main Execution ---
if __name__ == "__main__":
    answer_a = solve_part_a()
    answer_b = solve_part_b()
    answer_c = solve_part_c()
    
    final_answer = f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}]"
    
    print("\n--- Final Answer ---")
    print(f"<<<(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}]>>>")
