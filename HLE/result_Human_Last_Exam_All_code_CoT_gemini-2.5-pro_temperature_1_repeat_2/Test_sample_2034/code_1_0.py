import numpy as np

# Step 10: We have the identity c_1 * alpha_i - c_2 * (K*alpha)_i = alpha_i + (K^-1 * alpha)_i
# Step 11: For K close to I, K^-1 is approximately 2*I - K.
# Step 12: Substitute K^-1 into the identity.
# c_1 * alpha_i - c_2 * (K*alpha)_i = alpha_i + ((2*I - K)*alpha)_i
# c_1 * alpha_i - c_2 * (K*alpha)_i = alpha_i + 2*alpha_i - (K*alpha)_i
# c_1 * alpha_i - c_2 * (K*alpha)_i = 3*alpha_i - (K*alpha)_i
# This equation must hold for arbitrary alpha and K.
# By comparing the coefficients of alpha_i and (K*alpha)_i on both sides:
# Coefficient of alpha_i: c_1 = 3
# Coefficient of (K*alpha)_i: -c_2 = -1  => c_2 = 1

c1 = 3
c2 = 1

# Final equation form:
# -(K * alpha_D_minus_i)_i <= (1 + c1*beta)*alpha_D_i - (1 + c2*beta)*(K*alpha_D)_i + o(beta)
# Printing the final equation with the determined coefficients.
# Note: Python cannot print the exact mathematical notation, so we describe it.

print("The determined coefficients are:")
print(f"c1 = {c1}")
print(f"c2 = {c2}")
print("\nThe final equation is:")
print(f"-(K * a_D-i)_i <= (1 + {c1}*b)*a_D_i - (1 + {c2}*b)*(K*a_D)_i + o(b)")
print("where 'a' represents alpha and 'b' represents beta.")
