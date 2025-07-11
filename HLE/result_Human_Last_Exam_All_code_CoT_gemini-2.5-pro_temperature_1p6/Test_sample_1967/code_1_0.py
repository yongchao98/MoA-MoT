# The user can change this value for n.
n = 4

# The number of different sets of destabilizers for the n-qubit stabilizer
# set {Z_1, ..., Z_n} is given by the formula: 2**((n**2 + 5*n) / 2).

# Let's break down the calculation for the given n.

# The final formula is of the form: base**exponent.
# The numbers in the final equation are:
base = 2
coeff_n_squared = 1
coeff_n = 5
divisor = 2

# Calculate the exponent
# Using integer division // since the numerator is always even.
exponent_numerator = coeff_n_squared * n**2 + coeff_n * n
exponent = exponent_numerator // divisor

# Calculate the final result
total_sets = base**exponent

print(f"For a system with n = {n} qubits:")
print(f"The number of different destabilizer sets is calculated with the final equation:")
print(f"Result = {base}**(({coeff_n_squared}*n**2 + {coeff_n}*n) / {divisor})")
print("\nCalculation steps:")
print(f"Result = {base}**(({coeff_n_squared}*{n*n} + {coeff_n}*{n}) / {divisor})")
print(f"       = {base}**({exponent_numerator} / {divisor})")
print(f"       = {base}**{exponent}")
print(f"       = {total_sets}")