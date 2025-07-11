import collections

# Define the epitopes and their given IDs
epitopes = {
    "E1": "TYQRTRALV",
    "E2": "TFQRTRALV",
    "E3": "TFQRTRALK",
    "E4": "TYQRMFALV",
    "E5": "TYIPMRALV"
}

# The ranking is determined by analyzing the peptide sequences against the known
# binding motif of the H-2Kd MHC allele, as explained in the reasoning above.
# The primary factors are the anchor residues at positions P2 and P9.
#
# 1. E1 (TYQRTRALV): Reference high-affinity peptide. Optimal anchors (P2=Y, P9=V). HIGHEST affinity.
# 2. E5 (TYIPMRALV): Optimal anchors. Internal changes (Q->I, R->P, T->M). Less optimal than E1, but better than E4. 2nd HIGHEST affinity.
# 3. E4 (TYQRMFALV): Optimal anchors. Internal changes (T->M, R->F). R->F at P6 is unfavorable. 3rd HIGHEST affinity.
# 4. E2 (TFQRTRALV): Suboptimal P2 anchor (F instead of Y). Affinity is significantly reduced. 4th HIGHEST affinity.
# 5. E3 (TFQRTRALK): Catastrophic P9 anchor (K instead of V). Binding is nearly abrogated. LOWEST affinity.
#
# Final ranked order from highest to lowest affinity: E1, E5, E4, E2, E3

final_ranking = ["E1", "E5", "E4", "E2", "E3"]

# Print the final ranking
print("The epitopes ranked from highest to lowest expected amount complexed with H2-Kd are:")
# Using "equation" in the context of a ranked list, we output the final order.
print(", ".join(final_ranking))