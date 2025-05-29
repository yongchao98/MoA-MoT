def is_consecutive(seq):
    """Check if the given sequence of letters is consecutive."""
    return ord(seq[1]) == ord(seq[0]) + 1 and ord(seq[2]) == ord(seq[1]) + 1

def determine_result(cards, new_card):
    """Determine the result of adding a new card to the current hand."""
    # Check for "Peng"
    card_count = {}
    for card in cards:
        card_count[card] = card_count.get(card, 0) + 1
    if card_count.get(new_card, 0) >= 2:
        return 1  # Peng

    # Check for "Chi"
    cards.append(new_card)
    cards.sort()
    for i in range(len(cards) - 2):
        if is_consecutive(cards[i:i+3]):
            return 2  # Chi

    return 0  # Pass

# Initial cards
cards = list("FSNGLAWBLONTQ")

# Rounds
rounds = [
    ('H', 'A'),
    ('Y', 'B'),
    ('K', 'Q'),
    ('L', 'S'),
    ('N', None)  # No card to remove in the final round
]

# Process each round
for new_card, remove_card in rounds:
    if remove_card:
        cards.remove(remove_card)
    result = determine_result(cards, new_card)
    cards.append(new_card)
    cards.remove(new_card)  # Remove the added card to simulate the next round

# The result of the final round
print(result)