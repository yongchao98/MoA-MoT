import math

# For part (a), the reasoning leads to a "Yes".
answer_a = "Yes"

# For part (b), we calculate the number of similarity classes for n=3, with
# eigenvalues chosen from 3 distinct options. This is a combination with
# repetition problem.
n = 3  # The size of the matrix, which is the size of the multiset.
k = 3  # The number of distinct items to choose from.

# The formula for combinations with repetition is C(n + k - 1, n).
# Here, we calculate C(3 + 3 - 1, 3) = C(5, 3).
answer_b = math.comb(n + k - 1, n)

# For part (c), the reasoning shows that the growth is polynomial in n,
# not exponential.
answer_c = "No"

# We format the final answer as requested.
final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}."

print(final_answer)