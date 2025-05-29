# Initial cards
cards = list("WGLRQQYUWVLXZ")

# Define the rounds
rounds = [
    ('Q', 'W'),
    ('Q', 'V'),
    ('Z', 'W'),
    ('Q', 'X'),
    ('Q', None)  # No card to remove in the last round
]

# Function to determine the result of a round
def determine_result(cards, new_card):
    # Check for "Peng"
    if cards.count(new_card) >= 2:
        return 1  # Peng
    # Check for "Chi"
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return 2  # Chi
    return 0  # Pass

# Simulate the rounds
for new_card, remove_card in rounds:
    cards.append(new_card)
    if remove_card:
        cards.remove(remove_card)
    result = determine_result(cards, new_card)

# Print the result of the final round
print(result)