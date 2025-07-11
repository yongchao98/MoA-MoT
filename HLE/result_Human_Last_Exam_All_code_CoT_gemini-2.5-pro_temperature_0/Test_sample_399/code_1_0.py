from fractions import Fraction

# This script calculates the value of x(1) based on the derived particular solution
# x(t) = 1 / (t - t^2/2 - 1/8).

# We need to evaluate this at t=1.
t = Fraction(1)

# The equation for x(1) is 1 / (1 - 1^2/2 - 1/8)
val_t = t
val_t_squared_div_2 = t**2 / 2
constant_term = Fraction(1, 8)

# Calculate the denominator
denominator = val_t - val_t_squared_div_2 - constant_term

# Calculate the final result
result = 1 / denominator

print("The particular solution is x(t) = 1 / (t - t^2/2 - 1/8).")
print("To find x(1), we substitute t=1 into the equation.")
print(f"The final equation is: x(1) = 1 / ({val_t} - {val_t}**2/2 - {constant_term})")
print("\n--- Calculation Steps ---")
print(f"Denominator = {val_t} - {val_t_squared_div_2} - {constant_term}")
print(f"Denominator = {val_t - val_t_squared_div_2} - {constant_term}")
print(f"Denominator = {denominator}")
print(f"x(1) = 1 / {denominator}")
print(f"\nThe final value of x(1) is: {result}")