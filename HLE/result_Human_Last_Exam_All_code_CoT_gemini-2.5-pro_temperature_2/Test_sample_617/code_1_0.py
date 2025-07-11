#
# This script determines the color of a concentrated Pseudomonas aeruginosa sample
# based on the pigments it produces.
#

def find_bacterial_color():
    """
    Explains the color of a Pseudomonas aeruginosa pellet.
    """
    # Step 1: Define the primary pigments and their individual colors.
    # Pseudomonas aeruginosa is known for producing these characteristic pigments.
    pigment_1_name = "Pyocyanin"
    pigment_1_color = "Blue"

    pigment_2_name = "Pyoverdine"
    pigment_2_color = "Greenish-yellow"

    # Step 2: The final color of the dense bacterial sample is the result of
    # mixing these primary pigments.
    final_color = "Blue-green"
    final_choice = "C"

    # Step 3: Print the reasoning in the form of a conceptual equation.
    # Note: This is a qualitative description, not a mathematical equation.
    print("Determining the color of a concentrated Pseudomonas aeruginosa sample:")
    print(f"1. The bacterium produces a blue pigment called {pigment_1_name}.")
    print(f"2. It also produces a greenish-yellow pigment called {pigment_2_name}.")
    print("\nConceptual Color Equation:")
    print(f"Pigment Color '{pigment_1_color}' + Pigment Color '{pigment_2_color}' => Resulting Sample Color '{final_color}'")

    print(f"\nConclusion: The dense, washed pellet of Pseudomonas aeruginosa is {final_color}.")
    print(f"This corresponds to answer choice {final_choice}.")

# Execute the function to print the result.
find_bacterial_color()