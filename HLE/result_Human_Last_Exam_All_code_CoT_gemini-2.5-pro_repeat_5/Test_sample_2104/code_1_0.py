import math

# Step 1: Identify n1 and n2 based on the problem's conditions.
# As explained in the reasoning, the two smallest positive integers n
# satisfying the conditions are 6 and 8.
n1 = 6
n2 = 8

# Step 2: Calculate the parameter alpha.
alpha_numerator = n1 - 1
alpha_denominator = n2
alpha = alpha_numerator / alpha_denominator

# Step 3: Define and evaluate the function T(alpha).
# T(alpha) is identified as the Beta function B(alpha, 1-alpha),
# which is equal to pi / sin(pi * alpha).
# We need to calculate T(5/8) = pi / sin(5*pi/8).

# We will now compute this value.
# The numbers in the final equation are pi, 5, and 8.
val_pi = math.pi
val_num = 5
val_den = 8

# The calculation is pi / sin(5*pi/8)
result = val_pi / math.sin(val_num * val_pi / val_den)

# Print the values and the final equation as requested.
print(f"The 1st smallest integer is n1 = {n1}")
print(f"The 2nd smallest integer is n2 = {n2}")
print(f"The value alpha is ({n1} - 1) / {n2} = {alpha_numerator}/{alpha_denominator}")
print(f"We are asked to compute T({alpha_numerator}/{alpha_denominator})")
print(f"The function T(alpha) is interpreted as pi / sin(pi * alpha).")
print(f"The final equation is T({val_num}/{val_den}) = pi / sin({val_num}*pi/{val_den})")
print(f"The numerical result is: {result}")
