import matplotlib.pyplot as plt
import numpy as np

def plot_leucine_zipper_conservation():
    """
    This function generates and describes the conservation pattern
    of a canonical Leucine Zipper motif.
    """

    print("Analyzing the protein domain based on its conservation pattern:")
    print("The image displays a conservation plot with a clear periodic pattern. Highly conserved residues (red) appear at regular intervals.")
    print("This pattern is the hallmark of a coiled-coil structure, which is formed by alpha-helices that have a repeating sequence of seven amino acids called a heptad repeat (abcdefg).")
    print("\nThe Leucine Zipper is a classic example of a coiled-coil domain. Its key feature is the conservation of residues at specific positions in the heptad repeat.")
    print("Typically, positions 'a' and 'd' are conserved hydrophobic residues. This creates a pattern where conserved residues are separated by 2 and 3 other residues, respectively (a-x-x-d-x-x-x).")
    
    # Define the pattern for a canonical leucine zipper
    # Heptad repeat: a, b, c, d, e, f, g
    # Conserved positions: 'a' and 'd'
    sequence_length = 38
    # 1 for conserved, 0 for variable
    pattern = []
    for i in range(sequence_length):
        pos_in_heptad = (i % 7) + 1
        # Check if the position is 'a' (1) or 'd' (4)
        if pos_in_heptad == 1 or pos_in_heptad == 4:
            pattern.append(1)
        else:
            pattern.append(0)

    print("\nA simplified model for this pattern over 38 residues would be (1=conserved, 0=variable):")
    print(' '.join(map(str, pattern)))
    print("\nThis pattern of repeating conserved residues matches the fundamental characteristic shown in the image.")
    print("Other domains like SH3, PDZ, or Zinc Fingers have different, non-periodic conservation patterns related to their specific folds and functions.")
    print("\nTherefore, the domain shown is a Leucine zipper motif.")

# Execute the function to print the analysis
plot_leucine_zipper_conservation()