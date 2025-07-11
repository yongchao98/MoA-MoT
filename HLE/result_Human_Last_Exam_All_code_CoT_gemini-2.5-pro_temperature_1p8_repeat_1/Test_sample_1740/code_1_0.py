# The let-7 family is a group of microRNAs known for their role as tumor suppressors.
# The number of members is based on the distinct mature miRNA sequences identified in humans.
# This information is cataloged in miRBase, the central miRNA database.

# List of the distinct members of the let-7 family in humans (Homo sapiens),
# including the closely related miR-98.
let7_family_members = [
    "let-7a",
    "let-7b",
    "let-7c",
    "let-7d",
    "let-7e",
    "let-7f",
    "let-7g",
    "let-7i",
    "miR-98"
]

# Calculate the total number of members in the family.
total_members = len(let7_family_members)

# Print the final result.
print("The following members of the let-7 family have been identified in humans:")
for member in let7_family_members:
    print(f"- {member}")

print(f"\nTotal number of members: {let7_family_members[0]} + {let7_family_members[1]} + ... gives {total_members}")
print(f"There are {total_members} identified members in the human let-7 family.")
