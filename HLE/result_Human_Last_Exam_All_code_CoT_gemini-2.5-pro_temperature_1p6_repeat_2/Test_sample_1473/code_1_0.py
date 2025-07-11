import math

# The value of the integral is given by the expression I = pi * ln(1 + sqrt(2))
# Let's break down the components of this equation and compute the final value.

pi_val = math.pi
one_val = 1
sqrt2_val = math.sqrt(2)

# The final equation is I = pi * ln(1 + sqrt(2))
# We will print each numerical part of this equation.
print("The final equation is: I = pi * log(1 + sqrt(2))")
print(f"pi = {pi_val}")
print(f"{one_val} = {one_val}")
print(f"sqrt(2) = {sqrt2_val}")

# Now, we calculate the final result
final_value = pi_val * math.log(one_val + sqrt2_val)

print(f"\nThe final value of the integral is: {final_value}")
