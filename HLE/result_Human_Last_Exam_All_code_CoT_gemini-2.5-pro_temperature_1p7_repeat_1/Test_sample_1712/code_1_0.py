import math

# The problem is to evaluate the sum S = sum_{n=2 to inf} (-1)^n * S_n / n
# where S_n is the n-th harmonic number.

# Through analytical derivation, we found the closed-form expression for the sum:
# S = 1 + (ln(2))^2 / 2 - pi^2 / 12

print("The final closed-form expression for the sum is: S = 1 + (ln(2))^2 / 2 - pi^2 / 12")

# As requested, here are the numbers in the final equation, represented as:
# S = A + (ln(B))^C / D - pi^E / F
print("\nNumbers in the final equation:")
print("A = 1")
print("B = 2")
print("C = 2")
print("D = 2")
print("E = 2")
print("F = 12")

# We can verify this result by comparing the value from the formula
# with a numerical approximation of the sum.

# 1. Calculation using the formula
ln_2 = math.log(2)
pi = math.pi
exact_value = 1 + (ln_2**2) / 2 - (pi**2) / 12

print("\n--- Verification ---")
print(f"The exact value from the formula is: {exact_value}")

# 2. Numerical summation
def harmonic_number(n):
    """Calculates the n-th harmonic number."""
    # This is inefficient for large n, but sufficient for demonstration.
    # For very large n, approximations like ln(n) + gamma should be used.
    h = 0.0
    for k in range(1, n + 1):
        h += 1/k
    return h

def numerical_sum(num_terms):
    """Calculates the sum for a given number of terms."""
    total = 0.0
    current_h = 0.0
    for n in range(1, num_terms + 1):
        current_h += 1/n # More efficient calculation of harmonic numbers
        if n >= 2:
            total += ((-1)**n) * current_h / n
    return total

# Using a large number of terms for better accuracy
# Note: The convergence of this series is slow.
n_terms = 50000 
sum_approximation = numerical_sum(n_terms)

print(f"Numerical approximation with {n_terms} terms: {sum_approximation}")
print(f"The difference is: {abs(exact_value - sum_approximation)}")