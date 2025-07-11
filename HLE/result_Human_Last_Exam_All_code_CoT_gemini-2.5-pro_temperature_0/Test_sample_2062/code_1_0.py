import math

# (a) Is it true that A and B are similar if and only if the multiplicities of each eigenvalue are identical?
# For diagonal matrices, this is true. Similarity is determined by the multiset of eigenvalues.
answer_a = "Yes"

# (b) For n = 3 and 3 distinct eigenvalues, how many similarity classes exist?
# This is a stars and bars problem: choosing n items from k categories with replacement.
# n = 3 (matrix size)
# k = 3 (number of distinct eigenvalues available)
# The number of classes is given by the binomial coefficient C(n + k - 1, k - 1).
n_b = 3
k_b = 3
num_classes = math.comb(n_b + k_b - 1, k_b - 1)

# Format the expression for part (b) to show the calculation, as requested.
expression_b = f"C({n_b}+{k_b}-1, {k_b}-1) = C({n_b + k_b - 1}, {k_b - 1}) = {num_classes}"


# (c) Does the number of similarity classes grow exponentially with n for fixed q?
# The number of classes is C(n+q-1, q-1), which is a polynomial in n of degree q-1.
# Polynomial growth is not exponential.
answer_c = "No"

# Print the final answer in the required format.
print(f"(a) {answer_a}; (b) [{expression_b}]; (c) {answer_c}")
