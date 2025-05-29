from collections import Counter

# Initial cards
initial_cards = "ISZTWVHUUXWKC"

# New card added in Round 1
new_card = "P"

# Count occurrences of each card in the initial hand
card_count = Counter(initial_cards)

# Add the new card to the hand
card_count[new_card] += 1

# Check for the "Peng" condition
result = 0  # Default to "Pass"
for card, count in card_count.items():
    if count >= 3 and card == new_card:
        result = 1  # "Peng" condition met
        break

print(result)