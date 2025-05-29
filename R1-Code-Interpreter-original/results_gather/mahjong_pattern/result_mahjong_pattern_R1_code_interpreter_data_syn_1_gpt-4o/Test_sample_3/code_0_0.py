# Initial set of cards
cards = list("ISZTWVHUUXWKC")

# New card to be added
new_card = "P"

# Add the new card
cards.append(new_card)

# Check for "Peng" condition
result = 0  # Default to "Pass"
for card in set(cards):
    if cards.count(card) >= 3 and card == new_card:
        result = 1  # "Peng" condition met
        break

# Output the result
print(result)