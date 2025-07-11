import numpy as np
from math import factorial

# Plan: Demonstrate that if all moments M_k are zero, the Fourier transform must be zero.
#
# The Fourier transform F(ξ) can be expanded as a series involving the moments M_k:
# F(ξ) = Σ_{k=0 to ∞} [ c_k * M_k ], where the coefficients c_k = ((-2πiξ)^k / k!).
# The problem states that all moments M_k = ∫ x^k f(x) dx are 0.
# We will construct the equation for F(ξ) under this condition.

# Choose an arbitrary point ξ for demonstration purposes.
xi = 0.2

# The number of terms to display in our example equation.
num_terms = 5

print(f"The Fourier transform F(ξ) at ξ={xi} is given by the sum:")
print(f"F({xi}) = M_0 * c_0 + M_1 * c_1 + M_2 * c_2 + ...")
print("We are given that M_k = 0 for all k.\n")

# We build the equation term by term. Each term is zero because M_k is zero.
equation_terms = []
for k in range(num_terms):
    # The coefficient depends on xi and k.
    # We calculate it, but it will be multiplied by a moment M_k=0.
    coefficient = ((-2j * np.pi * xi)**k) / factorial(k)
    moment = 0
    term_value = coefficient * moment
    equation_terms.append(f"{term_value.real:.1f}")

final_result = 0

# Print the final equation with each number.
print("Final Equation:")
# The equation shows each term is 0, leading to a sum of 0.
equation_string = " + ".join(equation_terms)
print(f"F({xi}) = {equation_string} + ... = {final_result}")

print("\nSince ξ was chosen arbitrarily, F(ξ) must be 0 for all ξ.")
print("The injectivity of the Fourier transform on Schwartz functions implies that if F(ξ) = 0, then f(x) = 0.")