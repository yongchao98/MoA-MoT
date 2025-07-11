# This script classifies three graphene nanoribbon band structures
# based on their edge type, width, and electronic properties (metallic/semiconducting).

# --- Analysis of Plot 1 ---
# Edge: 'A' (Armchair) because of the symmetric band structure with a gap at k=0.
# Width: 8. There are 8 bands above E=0 and 8 below, so 2N=16 -> N=8.
# Band: '1' (Semiconducting) because of the clear energy gap at the Fermi level (E=0).
edge1 = "A"
width1 = 8
band1 = 1
classification1 = f"{edge1}{width1}{band1}"

# --- Analysis of Plot 2 ---
# Edge: 'Z' (Zigzag) due to the characteristic flat edge state at the Fermi level.
# Width: 6. There are 6 bands at or above E=0, so 2N=12 -> N=6.
# Band: '0' (Metallic) because there is no band gap.
edge2 = "Z"
width2 = 6
band2 = 0
classification2 = f"{edge2}{width2}{band2}"

# --- Analysis of Plot 3 ---
# Edge: 'A' (Armchair) because of the symmetric band structure.
# Width: 5. There are 5 bands at or above E=0, so 2N=10 -> N=5.
# Band: '0' (Metallic) because the valence and conduction bands touch at E=0.
edge3 = "A"
width3 = 5
band3 = 0
classification3 = f"{edge3}{width3}{band3}"

# Concatenate all classifications into a single string
final_classification = classification1 + classification2 + classification3

# Print the final result
print(final_classification)