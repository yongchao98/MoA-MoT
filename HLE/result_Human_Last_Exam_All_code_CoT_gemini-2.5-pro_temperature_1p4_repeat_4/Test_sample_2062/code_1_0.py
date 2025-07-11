import math
import numpy as np

# This script analyzes the similarity classes of diagonal matrices.

# --- Part (a) Analysis ---
# Question: Are two diagonal matrices A and B similar if and only if
# the multiplicities of each eigenvalue are identical?
#
# Explanation: Yes. The eigenvalues of a diagonal matrix are its diagonal
# entries. Two matrices are similar if and only if they have the same
# Jordan Normal Form. For a diagonal matrix, its Jordan form is a diagonal
# matrix with the eigenvalues sorted. Therefore, two diagonal matrices are
# similar if and only if their diagonal entries are permutations of each other,
# which is equivalent to having the same multiset of eigenvalues.

print("--- Part (a) Demonstration ---")
# Example matrices
A_diag = [1, 5, 5, 9]
B_diag = [5, 9, 1, 5] # A permutation of A's diagonal
C_diag = [1, 1, 5, 9] # Different multiplicities

# The eigenvalues are the diagonal elements. We check for similarity by
# comparing the sorted lists of eigenvalues.
eig_A = np.sort(A_diag)
eig_B = np.sort(B_diag)
eig_C = np.sort(C_diag)

print(f"Eigenvalues of A (sorted): {eig_A}")
print(f"Eigenvalues of B (sorted): {eig_B}")
print(f"Eigenvalues of C (sorted): {eig_C}")

# np.array_equal checks if two arrays have the same shape and elements.
are_A_B_similar = np.array_equal(eig_A, eig_B)
are_A_C_similar = np.array_equal(eig_A, eig_C)

print(f"\nA and B have the same multiset of eigenvalues. Are they similar? {are_A_B_similar}")
print(f"A and C have different multisets of eigenvalues. Are they similar? {are_A_C_similar}")
answer_a = "Yes"
print(f"Conclusion for (a): {answer_a}")


# --- Part (b) Analysis ---
# Question: For n=3 and distinct eigenvalues alpha, beta, gamma,
# how many similarity classes exist?
#
# Explanation: Since the eigenvalues are distinct, the multiset of eigenvalues
# for any such matrix is uniquely {alpha, beta, gamma}. By the logic of part (a),
# all diagonal matrices with these eigenvalues (e.g., diag(a,b,c), diag(b,a,c))
# belong to the same single similarity class.

print("\n--- Part (b) Analysis ---")
answer_b = 1
print(f"For n=3 and 3 distinct eigenvalues, there is only one possible multiset of eigenvalues.")
print(f"Therefore, there is only {answer_b} similarity class.")
print(f"Conclusion for (b): {answer_b}")


# --- Part (c) Analysis ---
# Question: Does the number of similarity classes for diagonal matrices in M_n(F)
# grow exponentially with n for a fixed q, where |F| = q?
#
# Explanation: The number of classes equals the number of multisets of size n
# chosen from q field elements. This is a "stars and bars" combinatorics problem,
# and the solution is C(n + q - 1, n). This formula is a polynomial in n of
# degree q-1. Polynomial growth is not exponential.

def combinations(n, k):
    """Helper function to calculate nCk (combinations)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def num_similarity_classes(n, q):
    """Calculates C(n + q - 1, n)."""
    return combinations(n + q - 1, n)

print("\n--- Part (c) Demonstration ---")
q_fixed = 5 # A fixed number of elements in the field F, e.g., F_5
print(f"We will check the growth rate for a fixed q = {q_fixed}.")
print("The number of classes is C(n + q - 1, n), a polynomial in n.")
print("We compare this to an exponential function, 2^n.")
print("\n" + f"{'n':>3} | {'Num Classes (Polynomial Growth)':>35} | {'Exponential Growth (2^n)':>25}")
print("-" * 70)
for n in range(1, 16):
    # C(n + 5 - 1, n) = C(n + 4, n)
    classes = num_similarity_classes(n, q_fixed)
    exponential_val = 2**n
    # The final equation C(n+4, n) = result is printed as requested
    print(f"{n:>3} | {f'C({n + q_fixed - 1}, {n}) = {classes}':>35} | {f'{exponential_val}':>25}")

print("\nThe number of classes grows polynomially (like n^4 here), much slower than 2^n.")
answer_c = "No"
print(f"Conclusion for (c): {answer_c}")

# --- Final Answer Summary ---
final_answer_string = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
print("\n" + "="*20 + " FINAL ANSWER " + "="*20)
print(final_answer_string)
print(f"<<<{final_answer_string}>>>")