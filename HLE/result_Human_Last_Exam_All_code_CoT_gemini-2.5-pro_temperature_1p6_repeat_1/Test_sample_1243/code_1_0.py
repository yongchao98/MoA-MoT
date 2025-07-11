import math

def solve_part_a():
    """
    Calculates the number of 2-isotropic vectors for defining a 2-neighbor of Z^16.
    This corresponds to counting non-zero vectors in (Z/2Z)^16 whose
    Hamming weight is a multiple of 4.
    """
    n = 16
    weights = [4, 8, 12, 16]
    terms = [math.comb(n, k) for k in weights]
    total = sum(terms)
    
    print("(a) The number of distinct 2-isotropic vectors is given by the sum:")
    equation = " + ".join([f"C({n}, {k})" for k in weights])
    values = " + ".join(map(str, terms))
    print(f"{equation} = {values} = {total}")
    return total

def solve_part_c():
    """
    Calculates the number of unimodular 2-neighbors of Z^12.
    This corresponds to counting non-zero vectors in (Z/2Z)^12 whose
    Hamming weight is a multiple of 4.
    """
    n = 12
    weights = [4, 8, 12]
    terms = [math.comb(n, k) for k in weights]
    total = sum(terms)

    print("\n(c) The number of unimodular 2-neighbors is given by the sum:")
    equation = " + ".join([f"C({n}, {k})" for k in weights])
    values = " + ".join(map(str, terms))
    print(f"{equation} = {values} = {total}")
    return total

# --- Main execution ---
# Calculate and print the result for part (a)
answer_a = solve_part_a()

# Explain the result for part (b)
answer_b = "yes"
print(f"\n(b) The answer is {answer_b.capitalize()}. It is possible for the neighbor to be even.")
print("This requires sum(x_i^2) to be a multiple of 18. A constructive example is the vector")
print("x = (1, 1, 1, 1, 1, 7, 6, 6) in Z^8, for which sum(x_i^2) = 126, which is 7 * 18.")

# Calculate and print the result for part (c)
answer_c = solve_part_c()

# Format the final answer string for submission
final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
# print(f"\nFinal Answer: {final_answer}")