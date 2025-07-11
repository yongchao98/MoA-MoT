import math

# Given binding affinities (macroscopic dissociation constants) in nM
Kd1 = 4.8
Kd2 = 11.2

# The relationship between Kd1, Kd2, and valency 'n' is:
# n = Kd2 / (Kd2 - 2 * Kd1)
# Let's calculate the terms in the equation.
numerator = Kd2
term_in_denominator = 2 * Kd1
denominator = Kd2 - term_in_denominator

# Calculate the valency 'n'
n_float = numerator / denominator
# Valency must be an integer
n = round(n_float)

print("To find the valency 'n', we solve the equation: n = Kd2 / (Kd2 - 2 * Kd1)\n")
print(f"Given values:")
print(f"Kd1 = {Kd1} nM")
print(f"Kd2 = {Kd2} nM\n")
print("Substituting the values into the equation:\n")
print(f"n = {numerator} / ({Kd2} - 2 * {Kd1})")
print(f"n = {numerator} / ({Kd2} - {term_in_denominator})")
print(f"n = {numerator} / {denominator}")
print(f"n = {n_float}\n")
print(f"Since valency must be an integer, the valency of the multimer is {n}.")
