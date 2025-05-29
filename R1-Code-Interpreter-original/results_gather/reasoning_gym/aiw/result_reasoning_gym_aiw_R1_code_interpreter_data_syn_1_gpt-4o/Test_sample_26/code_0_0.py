# Define the friends
linda_friends = {
    "Bob": {"Linda", "Tom", "Alice"},
    "Tom": {"Linda", "Bob", "Alice"},
    "Alice": {"Linda", "Bob", "Tom"}
}

# Count the number of female friends Bob has
bob_friends = linda_friends["Bob"]
female_friends_of_bob = {"Linda", "Alice"}  # Both are female

# Calculate the intersection to find how many female friends Bob has
female_friends_count = len(bob_friends.intersection(female_friends_of_bob))

print(female_friends_count)