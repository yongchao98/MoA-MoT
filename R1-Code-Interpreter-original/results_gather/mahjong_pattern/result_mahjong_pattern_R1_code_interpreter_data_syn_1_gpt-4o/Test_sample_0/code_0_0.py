# Initial cards
initial_cards = "QPDLJYCRXBSHN"

# New card to be added
new_card = "Y"

# Count occurrences of the new card in the initial cards
count_new_card = initial_cards.count(new_card)

# Determine the result
result = 1 if count_new_card >= 2 else 0

# Output the result
print(result)