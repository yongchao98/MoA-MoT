import math

def nCr_exact(n, r):
    """Calculates the binomial coefficient C(n, r) using math.comb for exact integers."""
    if r < 0 or r > n:
        return 0
    return math.comb(n, r)

# --- Part (a) Calculation ---
n_a = 16
k_values_a = [4, 8, 12, 16]
terms_a = [nCr_exact(n_a, k) for k in k_values_a]
result_a = sum(terms_a)

# --- Part (c) Calculation ---
n_c = 12
k_values_c = [4, 8, 12]
terms_c = [nCr_exact(n_c, k) for k in k_values_c]
result_c = sum(terms_c)

# --- Output ---
print("Part (a):")
print("The number of such vectors is the number of non-zero vectors in (Z/2Z)^16 whose Hamming weight is a multiple of 4.")
print("This is calculated by the sum of binomial coefficients C(16, k) for k = 4, 8, 12, 16.")
eq_a_str = " + ".join([f"C(16, {k})" for k in k_values_a])
val_a_str = " + ".join([str(t) for t in terms_a])
print(f"Number of vectors = {eq_a_str} = {val_a_str} = {result_a}")
print("-" * 20)

result_b_str = "no"
print("Part (b):")
print("For the neighbor lattice N_3(x) to be even, it requires that for all y in Z^8 with y.x = 0 (mod 3), y.y is even. This can be shown to lead to a contradiction, as it implies an all-odd vector must be an all-even vector. Therefore, it is not possible.")
print("-" * 20)

print("Part (c):")
print("The number of such neighbors is the number of non-zero vectors in (Z/2Z)^12 whose Hamming weight is a multiple of 4.")
print("This is calculated by the sum of binomial coefficients C(12, k) for k = 4, 8, 12.")
eq_c_str = " + ".join([f"C(12, {k})" for k in k_values_c])
val_c_str = " + ".join([str(t) for t in terms_c])
print(f"Number of neighbors = {eq_c_str} = {val_c_str} = {result_c}")
print("-" * 20)

final_answer_str = f"(a) {result_a}; (b) {result_b_str}; (c) {result_c}"
print(f"Final Answer: {final_answer_str}")

print(f"<<<{final_answer_str}>>>")