# This script classifies three graphene nanoribbon band structures based on visual analysis.

# --- Analysis of Ribbon 1 ---
# Edge: Symmetric around k=0 -> Armchair ('A')
# Width: 8 bands above E=0, 8 below -> Total 16 bands = 2*N -> N=8
# Band: Clear gap at E=0 -> Semiconducting ('1')
ribbon1_edge = 'A'
ribbon1_width = 8
ribbon1_band = 1
classification1 = f"{ribbon1_edge}{ribbon1_width}{ribbon1_band}"

# --- Analysis of Ribbon 2 ---
# Edge: Flat edge states near E=0 -> Zigzag ('Z')
# Width: 5 bands above E=0, 5 below -> Total 10 bands = 2*N -> N=5
# Band: Bands cross E=0 -> Metallic ('0')
ribbon2_edge = 'Z'
ribbon2_width = 5
ribbon2_band = 0
classification2 = f"{ribbon2_edge}{ribbon2_width}{ribbon2_band}"

# --- Analysis of Ribbon 3 ---
# Edge: Symmetric around k=0 -> Armchair ('A')
# Width: 6 bands above E=0, 6 below -> Total 12 bands = 2*N -> N=6
# Band: Bands touch at E=0 -> Metallic ('0')
ribbon3_edge = 'A'
ribbon3_width = 6
ribbon3_band = 0
classification3 = f"{ribbon3_edge}{ribbon3_width}{ribbon3_band}"

# Concatenate all classifications into the final string
final_classification = classification1 + classification2 + classification3

# Print the final result
print(final_classification)