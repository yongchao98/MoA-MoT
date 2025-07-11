# The let-7 family in humans is a group of microRNAs sharing a conserved seed sequence.
# The members are well-established in biological databases.
# This script lists the primary members and then prints the total count.

# List of human let-7 family members, including miR-98 which is part of the family.
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

# Calculate the number of members in the list
number_of_members = len(let7_family_members)

# Print the final answer
print(f"To date, {number_of_members} members of the let-7 family have been identified in humans.")