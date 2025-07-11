import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# The number of basis vectors is 6
n = 6

# The vectors are constructed from pairs of basis vectors
k = 2

# Calculate the number of pairs of indices
num_pairs = combinations(n, k)

# For each pair, there are two possible vectors (from the sign choice)
num_sign_choices = 2

# Calculate the total number of vectors
total_vectors = num_pairs * num_sign_choices

# Print the step-by-step calculation of the final equation
print(f"The number of ways to choose 2 distinct basis vectors from 6 is C(6, 2).")
print(f"C(6, 2) = (6 * 5) / (2 * 1) = {num_pairs}")
print(f"For each pair, there are {num_sign_choices} choices for the sign, leading to two distinct vectors.")
print(f"The total number of vectors is the product of these two numbers.")
print(f"{num_pairs} * {num_sign_choices} = {total_vectors}")
