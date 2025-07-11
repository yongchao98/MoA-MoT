# Define the composition based on a representative optimal value of 15 atomic percent (at%) Nickel.
# This means for every 100 total metal atoms, 15 are Ni and 85 are Ce.
ni_atoms = 15
ce_atoms = 85
total_atoms = ni_atoms + ce_atoms

# Calculate the simple atomic ratio of Ni to Ce
ni_ce_ratio = ni_atoms / ce_atoms

# Calculate the atomic ratio of Ni to the total metal content (Ni+Ce)
# This is another common way to express the ratio and corresponds to the atomic percent.
ni_total_ratio = ni_atoms / total_atoms

print("This script calculates the Ni/Ce ratio based on a representative optimal composition.")
print("A common optimal composition is around 15 atomic percent (at%) Ni.\n")
print("For this composition, out of 100 total metal atoms:")
# The prompt requires printing each number in the final equation.
print(f"Number of Ni atoms = {ni_atoms}")
print(f"Number of Ce atoms = {ce_atoms}\n")

print("The corresponding simple atomic ratio (Ni/Ce) is calculated as:")
print(f"{ni_atoms} / {ce_atoms} = {ni_ce_ratio:.4f}\n")

print("The atomic ratio of Ni to total metal content (Ni/(Ni+Ce)) is:")
print(f"{ni_atoms} / ({ni_atoms} + {ce_atoms}) = {ni_total_ratio:.2f}")
print("\nThis Ni/(Ni+Ce) ratio of 0.15 is often cited as being near the ideal value for maximizing catalytic performance.")
