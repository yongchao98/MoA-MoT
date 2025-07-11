import collections

# Proposing a reasonable set of partial charges for a new all-atom methanol model.
# This set must be neutral, treat the methyl hydrogens as equivalent, and
# be chemically sensible based on electronegativity.

# Choice D is selected as it meets all criteria.
charges = collections.OrderedDict({
    "Carbon": 0.1450,
    "Methyl H1": 0.0400,
    "Methyl H2": 0.0400,
    "Methyl H3": 0.0400,
    "Oxygen": -0.6830,
    "Hydroxyl hydrogen": 0.4180
})

print("Proposed reasonable partial charges for a new methanol model:")
for atom, charge in charges.items():
    # Using .ljust() for cleaner alignment
    print(f"{atom.ljust(18)}:\t{charge:+.4f}")

# Verify that the total charge of the molecule is zero and print the equation.
total_charge = sum(charges.values())

print("\nVerifying charge neutrality of the molecule:")
# Build the equation string with each number
charge_values = [f"({charge:+.4f})" for charge in charges.values()]
equation_str = " + ".join(charge_values)

# We use round() to handle potential floating point inaccuracies for the final sum
print(f"{equation_str} = {round(total_charge, 4)}")
print("\nThe sum of the partial charges is zero, confirming the molecule is neutral.")
