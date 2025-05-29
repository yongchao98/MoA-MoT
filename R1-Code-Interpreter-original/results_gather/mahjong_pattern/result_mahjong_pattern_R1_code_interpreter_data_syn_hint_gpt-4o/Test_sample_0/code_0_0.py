# Initial set of cards
initial_cards = "QPDLJYCRXBSHN"

# New card to be added
new_card = "Y"

# Count the occurrences of 'Y' in the initial set
count_Y = initial_cards.count(new_card)

# Determine the result based on the count
result = 1 if count_Y >= 2 else 0

# Print the result
print(result)