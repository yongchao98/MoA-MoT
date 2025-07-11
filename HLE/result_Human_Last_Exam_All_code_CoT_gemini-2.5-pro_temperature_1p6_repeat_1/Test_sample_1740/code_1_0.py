# The let-7 family is a group of microRNAs that are highly conserved across species.
# In humans, this family is composed of several members which share a similar seed sequence.
# The primary members are let-7a, let-7b, let-7c, let-7d, let-7e, let-7f, let-7g, let-7i, and miR-98.
# Some members, like let-7a and let-7f, are encoded by multiple genes but produce identical mature microRNAs.
# This script lists the unique mature members and calculates the total number.

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

# Calculate the total number of members
num_members = len(let7_family_members)

# Print the members
print("The members of the let-7 family in humans are:")
for member in let7_family_members:
    print(f"- {member}")

# Print the final count
print(f"\nTo date, {num_members} members of the let-7 family have been identified in humans.")
