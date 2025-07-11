# A list of the human let-7 family precursor miRNA genes (paralogs).
# This list is based on data from miRBase, the primary miRNA sequence database.
# The family includes let-7 paralogs and mir-98, which shares the same seed sequence family.
human_let7_family_genes = [
    "hsa-let-7a-1",
    "hsa-let-7a-2",
    "hsa-let-7a-3",
    "hsa-let-7b",
    "hsa-let-7c",
    "hsa-let-7d",
    "hsa-let-7e",
    "hsa-let-7f-1",
    "hsa-let-7f-2",
    "hsa-let-7g",
    "hsa-let-7i",
    "hsa-mir-98"
]

# Calculate the total number of members by getting the length of the list.
total_members = len(human_let7_family_genes)

print("The following are the known precursor genes for the let-7 family in humans:")

# Print each member that contributes to the final count.
# This represents the "equation" of all the individual members adding up to the total.
for i, member in enumerate(human_let7_family_genes):
    # This loop shows each component of our final count.
    # The first item is 1, the second is 2, and so on, up to the total.
    if i < len(human_let7_family_genes) - 1:
        print(f"1 ({member}) +")
    else:
        print(f"1 ({member})")

print(f"\n= {total_members}")
print(f"\nTo date, {total_members} members of the let-7 family have been identified in humans.")