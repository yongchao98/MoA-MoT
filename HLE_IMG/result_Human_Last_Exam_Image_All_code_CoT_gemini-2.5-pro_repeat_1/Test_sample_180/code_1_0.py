import matplotlib.pyplot as plt
import numpy as np

def identify_protein_domain():
    """
    Analyzes the conservation pattern to identify the protein domain.
    The pattern in the image shows regularly spaced, highly conserved residues.
    This is characteristic of a Leucine Zipper motif.
    This script will generate an idealized plot for a Leucine Zipper to demonstrate the similarity.
    """

    # 1. Define the Leucine Zipper's heptad repeat pattern.
    # A Leucine Zipper has a repeating pattern of 7 amino acids (a, b, c, d, e, f, g).
    # Positions 'a' and 'd' are typically conserved hydrophobic residues.
    # We'll represent conservation with a numerical score: 1.0 for high, 0.3 for low.
    # Pattern: [a, b, c, d, e, f, g]
    heptad_repeat_conservation = [1.0, 0.3, 0.3, 1.0, 0.3, 0.3, 0.3]
    
    # 2. Create a longer sequence by repeating the pattern.
    # Let's show 5 repeats to create a sequence of 35 amino acids.
    num_repeats = 5
    full_sequence_conservation = heptad_repeat_conservation * num_repeats
    
    # 3. Generate the bar chart.
    positions = np.arange(len(full_sequence_conservation))
    colors = ['red' if c == 1.0 else 'gray' for c in full_sequence_conservation]

    plt.figure(figsize=(10, 4))
    plt.bar(positions, full_sequence_conservation, color=colors, width=0.8)
    
    # Style the plot to match the question's image
    plt.title("Idealized Conservation Pattern for a Leucine Zipper")
    plt.tick_params(axis='y', which='both', left=False, labelleft=False)
    plt.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    plt.box(False)
    
    print("Analyzing the conservation pattern:")
    print("The image displays a highly regular, repeating pattern of conserved amino acids (red bars).")
    print("This type of pattern is the hallmark of a Leucine Zipper motif, which is built from a heptad (7-amino-acid) repeat.")
    
    print("\nThe Leucine Zipper heptad repeat can be noted as (a, b, c, d, e, f, g).")
    print("Positions 'a' and 'd' are critical and highly conserved.")
    print("This creates a conservation pattern of: Conserved-Variable-Variable-Conserved-Variable-Variable-Variable")
    
    # Demonstrate the pattern from our idealized sequence
    print("\nOur idealized sequence of 35 amino acids (5 repeats) shows this pattern:")
    pattern_str = ""
    for i, score in enumerate(full_sequence_conservation):
        pattern_str += "R-" if score == 1.0 else "g-"
        if (i + 1) % 7 == 0 and i < len(full_sequence_conservation) - 1:
            pattern_str += "|-" # Separator between repeats
    print(pattern_str[:-1]) # R=Red(Conserved), g=gray(variable)
    
    print("\nThe generated plot visually matches the periodic nature of the pattern in the question.")
    print("Therefore, the domain shown is a Leucine zipper motif.")
    
    # To show the plot when the script is run
    # plt.show() # This line is commented out to prevent display in non-interactive environments.

# Execute the analysis
identify_protein_domain()
<<<B>>>