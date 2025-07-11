# This script classifies three graphene nanoribbon band structures
# based on their visual features.

# Classification for Ribbon 1
# Edge: 'A' (Armchair), determined by the gap centered at k=0.
# Width: 8, determined by counting 8 bands above the Fermi level.
# Band: 1 (Semiconducting), determined by the presence of a band gap.
edge1 = 'A'
width1 = 8
band1 = 1

# Classification for Ribbon 2
# Edge: 'Z' (Zigzag), determined by flat bands and crossing away from k=0.
# Width: 7, determined by counting 7 bands above the Fermi level.
# Band: 0 (Metallic), determined by the absence of a band gap.
edge2 = 'Z'
width2 = 7
band2 = 0

# Classification for Ribbon 3
# Edge: 'A' (Armchair), determined by the bands meeting at k=0.
# Width: 9, determined by counting 9 bands above the Fermi level.
# Band: 0 (Metallic), determined by the absence of a band gap.
edge3 = 'A'
width3 = 9
band3 = 0

# The problem asks to output each number in the final equation.
# The code below builds the final string from these numbers and letters.
classification1 = f"{edge1}{width1}{band1}"
classification2 = f"{edge2}{width2}{band2}"
classification3 = f"{edge3}{width3}{band3}"

# Concatenate the classifications for the final result
final_string = classification1 + classification2 + classification3

print(final_string)