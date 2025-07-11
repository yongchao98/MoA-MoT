# Define the properties of the peptide synthesis
peptide_length = 100
# The number of coupling reactions is one less than the length
num_couplings = peptide_length - 1

# A list of possible coupling efficiencies per step
step_efficiencies = [0.995, 0.99, 0.98, 0.97]

print(f"Calculating theoretical yield for a {peptide_length}aa peptide after {num_couplings} coupling steps:\n")

for efficiency in step_efficiencies:
    # The formula for total yield is (step_efficiency) ^ (number_of_couplings)
    total_yield = efficiency ** num_couplings

    # Print the result, showing each number in the equation
    print(
        f"With {efficiency:.1%} efficiency per step: "
        f"({efficiency}) ** {num_couplings} = {total_yield:.2%} final yield"
    )
