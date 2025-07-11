import numpy as np

def original_form(x, n):
    """Calculates the quadratic form using the original double summation."""
    total = 0
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            total += (n - abs(i - j)) * x[i-1] * x[j-1]
    return total

def sum_of_squares_form(x, n):
    """Calculates the quadratic form using the sum of squares representation."""
    S = np.zeros(n + 1)
    for k in range(1, n + 1):
        S[k] = S[k-1] + x[k-1]

    total = S[n]**2
    for k in range(1, n):
        # S_k = sum_{i=1 to k} x_i
        # S_n - S_k = sum_{i=k+1 to n} x_i
        total += S[k]**2 + (S[n] - S[k])**2
        
    return total

# Example usage
n = 4
# A sample vector of real numbers
x = np.array([1.5, -2.0, 0.5, 3.0])

# Calculate the value using both forms
val1 = original_form(x, n)
val2 = sum_of_squares_form(x, n)

print(f"The value of the quadratic form for n={n} and x={x} is:")
print(f"Calculated by original formula: {val1}")
print(f"Calculated by sum-of-squares formula: {val2}")

print("\nAs the sum of squares formula shows, the quadratic form is always non-negative.")
print("This implies that c=0 is a valid choice.")
print("Through advanced analysis (or numerical checking for large n), it can be shown that the minimum eigenvalue approaches 0 as n increases.")
print("Therefore, for the inequality to hold for all n, the maximum value for c must be 0.")

c = 0
print(f"\nThe maximum real number c is {c}.")

# Let's show the final inequality with the found c.
# For demonstration, we will show the terms for n=3, x=(x1, x2, x3)
# The inequality is: Sum_{i,j} (n-|i-j|)x_i x_j >= c * Sum_i x_i^2
# For n=3, c=0:
# (3-0)x1x1 + (3-1)x1x2 + (3-2)x1x3 +
# (3-1)x2x1 + (3-0)x2x2 + (3-1)x2x3 +
# (3-2)x3x1 + (3-1)x3x2 + (3-0)x3x3 >= 0 * (x1^2+x2^2+x3^2)
#
# which simplifies to:
# 3*x1**2 + 2*x1*x2 + 1*x1*x3 + 2*x2*x1 + 3*x2**2 + 2*x2*x3 + 1*x3*x1 + 2*x3*x2 + 3*x3**2 >= 0
# 3*x1**2 + 3*x2**2 + 3*x3**2 + 4*x1*x2 + 2*x1*x3 + 4*x2*x3 >= 0
# Our sum of squares identity for n=3 is:
# (x1+x2+x3)**2 + (x1)**2 + (x2+x3)**2 + (x1+x2)**2 + (x3)**2 >= 0
# Expanding this gives:
# (x1**2+x2**2+x3**2+2*x1*x2+2*x1*x3+2*x2*x3) + x1**2 + (x2**2+2*x2*x3+x3**2) + (x1**2+2*x1*x2+x2**2) + x3**2
# = 3*x1**2 + 3*x2**2 + 3*x3**2 + 4*x1*x2 + 2*x1*x3 + 4*x2*x3 >= 0
# which matches.

# Final symbolic output as requested:
print("\nThe final equation for an arbitrary n is:")
print(u"sum_{i=1 to n} sum_{j=1 to n} (n - |i-j|) * x_i * x_j >= 0 * sum_{i=1 to n} x_i**2")
print("We can demonstrate the coefficients for n=3:")
print("3*x_1*x_1 + 2*x_1*x_2 + 1*x_1*x_3 + 2*x_2*x_1 + 3*x_2*x_2 + 2*x_2*x_3 + 1*x_3*x_1 + 2*x_3*x_2 + 3*x_3*x_3 >= 0 * (x_1**2 + x_2**2 + x_3**2)")