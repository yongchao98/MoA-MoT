# Initial cards
initial_cards = "PFVYBJBNLYUCT"

# New card added in Round 1
new_card = "T"

# Add the new card to the hand
current_hand = initial_cards + new_card

# Check for "Peng"
result = 0  # Default to "Pass"
for card in set(current_hand):
    if current_hand.count(card) >= 3 and card == new_card:
        result = 1  # "Peng"
        break

print(result)