# Define the optimal atomic percentage of Ni based on experimental literature findings.
# A value around 15% is often cited for maximizing the catalytically active Ni-CeO2 interface.
ni_atomic_percent = 15.0

# The remainder of the composition is Cerium (Ce).
ce_atomic_percent = 100.0 - ni_atomic_percent

# The molar ratio is the ratio of their atomic percentages.
# We can think of this as considering a sample of 100 atoms.
ni_moles = ni_atomic_percent
ce_moles = ce_atomic_percent

# Calculate the final Ni/Ce ratio.
ideal_ratio = ni_moles / ce_moles

print(f"Based on an optimal composition of {ni_atomic_percent}% Ni and {ce_atomic_percent}% Ce (atomic percent):")
print("\nThe calculation for the ideal Ni/Ce molar ratio is as follows:")
print(f"Ni moles = {ni_moles}")
print(f"Ce moles = {ce_moles}")
print(f"Ratio Equation: {ni_moles} / {ce_moles} = {ideal_ratio:.3f}")

print(f"\nAn ideal Ni/Ce molar ratio is approximately {ideal_ratio:.3f}.")