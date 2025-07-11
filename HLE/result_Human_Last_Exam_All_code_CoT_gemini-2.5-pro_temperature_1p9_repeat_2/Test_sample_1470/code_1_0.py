import math

# Propose reasonable partial charges for the atoms of methanol (CH3OH).
# Based on chemical principles and comparison with established force fields,
# the charges from choice D are the most reasonable.

# Atom names and their corresponding proposed charges.
atoms = [
    ("Carbon", 0.1450),
    ("Methyl H1", 0.0400),
    ("Methyl H2", 0.0400),
    ("Methyl H3", 0.0400),
    ("Oxygen", -0.6830),
    ("Hydroxyl hydrogen", 0.4180)
]

print("Proposed Partial Charges for Methanol (CH3OH):")
print("-" * 45)

# Initialize variables for the verification step
total_charge = 0
charge_list = []

for atom, charge in atoms:
    print(f"{atom:<20} {charge:>8.4f}")
    total_charge += charge
    charge_list.append(charge)

print("-" * 45)
print("Verification of Net Charge:")

# Build the equation string
equation_parts = []
for charge in charge_list:
    if charge >= 0:
        equation_parts.append(f"{charge:.4f}")
    else:
        # Enclose negative numbers in parentheses for clarity
        equation_parts.append(f"({charge:.4f})")

equation_str = " + ".join(equation_parts)

# Use math.isclose() for robust floating-point comparison to zero
if math.isclose(total_charge, 0.0, abs_tol=1e-9):
    result = 0.0
else:
    result = total_charge
    
print(f"Equation: {equation_str} = {result:.4f}")

if math.isclose(result, 0.0, abs_tol=1e-9):
    print("The total charge is zero. The charge set is valid.")
else:
    print("Warning: The total charge is not zero!")
