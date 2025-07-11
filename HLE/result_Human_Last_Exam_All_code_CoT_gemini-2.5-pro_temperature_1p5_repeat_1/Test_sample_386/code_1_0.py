def print_culture_conditions():
    """
    This function defines and prints the key numerical parameters
    for the successful culture medium described in the correct answer.
    The correct procedure involves culturing stromal cells in a medium
    containing 10% FBS and 1% antibiotic.
    """

    # Percentage of Fetal Bovine Serum (FBS) in the medium
    fbs_percentage = 10

    # Percentage of antibiotic in the medium
    antibiotic_percentage = 1

    # Print the final equation representing the key supplements
    print("For successful proliferation, the medium was supplemented as follows:")
    print(f"Final Medium = Basal Medium + {fbs_percentage}% FBS + {antibiotic_percentage}% Antibiotic")

# Execute the function to display the output
print_culture_conditions()