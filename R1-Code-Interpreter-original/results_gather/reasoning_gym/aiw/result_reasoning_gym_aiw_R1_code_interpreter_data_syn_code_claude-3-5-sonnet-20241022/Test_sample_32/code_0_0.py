# Define the groups
male_friends = 4  # Total male friends of Alice
female_friends = 5  # Total female friends of Alice

# Calculate Charles's female friends
# Charles is friends with:
# - All of Alice's female friends (5)
# - Alice herself (1)
charles_female_friends = female_friends + 1  # Adding Alice

print(f"Charles has {charles_female_friends} female friends")