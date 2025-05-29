# Define the people
margaret = "Margaret"
james = "James"
other_male = "Male2"
female_friend = "Female1"

# Create a dictionary to store friendships
friends = {
    margaret: [james, other_male, female_friend],
    james: [margaret, other_male, female_friend],
    other_male: [margaret, james, female_friend],
    female_friend: [margaret, james, other_male]
}

# Count James's female friends
james_female_friends = sum(1 for friend in friends[james] if friend in [margaret, female_friend])
print(f"James has {james_female_friends} female friends")