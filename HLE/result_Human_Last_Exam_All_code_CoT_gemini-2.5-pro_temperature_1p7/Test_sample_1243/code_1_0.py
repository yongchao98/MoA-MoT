import math

# Part (a)
# We need to compute the number of non-zero vectors in F_2^16 with weight divisible by 4.
# This is C(16, 4) + C(16, 8) + C(16, 12) + C(16, 16).
n_a = 16
weights_a = [4, 8, 12, 16]
terms_a = [math.comb(n_a, k) for k in weights_a]
result_a = sum(terms_a)
equation_a_str = " + ".join(map(str, terms_a))

# Part (b)
# This is a theoretical question.
answer_b = "no"

# Part (c)
# We need to compute the number of non-zero vectors in F_2^12 with weight divisible by 4.
# This is C(12, 4) + C(12, 8) + C(12, 12).
n_c = 12
weights_c = [4, 8, 12]
terms_c = [math.comb(n_c, k) for k in weights_c]
result_c = sum(terms_c)
equation_c_str = " + ".join(map(str, terms_c))

# Print the detailed breakdown of calculations and the final answer.
print(f"For part (a), the number of vectors is the sum of C(16, k) for k in {{4, 8, 12, 16}}.")
print(f"The calculation is: {equation_a_str} = {result_a}")
print("-" * 20)
print(f"For part (b), the answer is {answer_b}.")
print("Explanation: The lattice Z^8 is odd. For d=3 (>2), its neighbors must also be odd and cannot be even.")
print("-" * 20)
print(f"For part (c), the number of neighbors is the sum of C(12, k) for k in {{4, 8, 12}}.")
print(f"The calculation is: {equation_c_str} = {result_c}")
print("-" * 20)
print(f"Final combined answer:")
print(f"(a) {result_a}; (b) {answer_b}; (c) {result_c}.")
