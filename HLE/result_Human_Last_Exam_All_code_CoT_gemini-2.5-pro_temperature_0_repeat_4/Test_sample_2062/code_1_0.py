import math

# (a) Is it true that A and B are similar if and only if the multiplicities of each eigenvalue are identical?
# For diagonal matrices, similarity is equivalent to having the same multiset of diagonal entries.
# This is the same as having identical multiplicities for each eigenvalue.
answer_a = "Yes"

# (b) For n = 3 and eigenvalues alpha, beta, gamma in F, how many similarity classes exist if alpha, beta, and gamma are distinct?
# A similarity class is defined by the multiset of eigenvalues. We need to find the number of
# multisets of size k=3 from a set of n=3 distinct elements {alpha, beta, gamma}.
# The formula is the multiset coefficient C(n+k-1, k).
n = 3
k = 3
# Calculate C(3+3-1, 3) = C(5, 3)
num_classes = math.comb(n + k - 1, k)
answer_b = num_classes

# (c) Does the number of similarity classes for diagonal matrices in M_n(F) grow exponentially with n for fixed q?
# The number of similarity classes is the number of multisets of size n from the q elements of the field F.
# This is given by C(n+q-1, n), which is a polynomial in n of degree q-1.
# Polynomial growth is not exponential.
answer_c = "No"

# Print the final answer in the specified format.
print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")