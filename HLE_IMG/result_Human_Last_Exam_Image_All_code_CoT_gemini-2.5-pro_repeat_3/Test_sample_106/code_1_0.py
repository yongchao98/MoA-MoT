# This script determines the identity of substituents after a backbone rearrangement.
# The analysis is based on the well-established Wagner-Meerwein rearrangement mechanism
# for pentacyclic triterpenoids.

# Step 1: Define a dictionary to hold the results of the chemical analysis.
# The keys are the numbered positions, and the values are the corresponding substituents.
substituents = {
    1: 'CH3',
    2: 'H',
    3: 'CH3',
    4: 'CH3',
    5: 'CH3'
}

# Step 2: Create a list to store the formatted strings for the output.
output_strings = []

# Step 3: Iterate through the dictionary in numerical order of the positions.
# The `sorted()` function ensures the output is ordered from 1 to 5.
for position in sorted(substituents.keys()):
    substituent = substituents[position]
    # Format each result as "position = substituent"
    output_strings.append(f"{position} = {substituent}")

# Step 4: Join the individual strings with ", " and print the final result.
print(", ".join(output_strings))