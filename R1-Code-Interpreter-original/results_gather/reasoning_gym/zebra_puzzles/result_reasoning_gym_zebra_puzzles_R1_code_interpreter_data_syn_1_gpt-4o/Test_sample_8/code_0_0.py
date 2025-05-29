# Initialize the houses with empty values
houses = [
    {"name": None, "smoothie": None, "lunch": None, "drink": None},
    {"name": None, "smoothie": None, "lunch": None, "drink": None},
    {"name": None, "smoothie": None, "lunch": None, "drink": None},
    {"name": None, "smoothie": None, "lunch": None, "drink": None}
]

# Apply the clues
houses[2]["name"] = "arnold"
houses[2]["lunch"] = "stir fry"
houses[2]["drink"] = "milk"
houses[2]["smoothie"] = "butterscotch"

houses[1]["name"] = "alice"

houses[0]["smoothie"] = "desert"

# Deduce the remaining information
# Since Arnold is in the third house, Alice is in the second house, and the Desert smoothie lover is in the first house,
# Bob must be in the fourth house, and Carol must be in the first house.

houses[3]["name"] = "bob"
houses[0]["name"] = "carol"

# From clue 4, Bob and the Dragonfruit smoothie lover are next to each other.
# Since Bob is in the fourth house, the Dragonfruit smoothie lover must be in the third house.

houses[3]["smoothie"] = "dragonfruit"

# From clue 5, the pizza lover drinks water. Since Arnold loves stir fry, Alice must love pizza and drink water.

houses[1]["lunch"] = "pizza"
houses[1]["drink"] = "water"

# From clue 8, the grilled cheese lover is directly left of the boba tea drinker.
# Since Alice is in the second house and drinks water, Carol must love grilled cheese and be in the first house.
# Therefore, Bob must drink boba tea.

houses[0]["lunch"] = "grilled cheese"
houses[3]["drink"] = "boba tea"

# The only remaining lunch is soup, which must be Bob's favorite.
houses[3]["lunch"] = "soup"

# The only remaining smoothie is darkness, which must be Alice's favorite.
houses[1]["smoothie"] = "darkness"

# Print the name of the person in the first house
print(houses[0]["name"])