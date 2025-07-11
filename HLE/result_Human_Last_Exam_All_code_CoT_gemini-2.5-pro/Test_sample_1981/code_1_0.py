import math

# Step 1: Define parameters from the problem description
# p=6, q=4, A=3, B=2A=6, T=ln(10), alpha=0
# These are used to derive the relation between y0 and x0 and the value of X0.
# The derivation shows that X0**(15/2) = 10**117.

# Step 2: Calculate X0
# X0**(15/2) = 10**117
# X0 = (10**117)**(2/15) = 10**(117 * 2 / 15) = 10**(234/15) = 10**15.6
X0 = 10**15.6

# Step 3: Calculate the terms needed for the final expression
X0_sq = X0**2
coeff1 = 10**30
coeff2 = -10**30
const = 10

# Step 4: Calculate the final result
result = coeff1 * X0_sq + coeff2 * X0 + const

# Step 5: Print the final equation with the computed values
# Using scientific notation for clarity
print(f"The value of X0 is approximately {X0:.4e}")
print(f"The value of X0^2 is approximately {X0_sq:.4e}")
print("\nThe final equation is:")
print(f"{coeff1:.1e} * {X0_sq:.4e} + ({coeff2:.1e}) * {X0:.4e} + {const}")
print(f"= {result:.4e}")

# The final answer is the numerical result of the expression.
final_answer = result
# The output format should be <<<answer>>>
# print(f"<<<{final_answer:.4e}>>>")