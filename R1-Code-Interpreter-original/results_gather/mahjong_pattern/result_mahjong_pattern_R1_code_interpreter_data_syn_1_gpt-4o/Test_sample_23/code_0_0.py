def is_consecutive(seq):
    """Check if the sequence of characters is consecutive."""
    return ord(seq[1]) == ord(seq[0]) + 1 and ord(seq[2]) == ord(seq[1]) + 1

def determine_result(cards, new_card):
    """Determine the result after adding a new card."""
    # Check for Peng
    for card in set(cards):
        if cards.count(card) == 2 and new_card == card:
            return 1  # Peng

    # Check for Chi
    cards.append(new_card)
    cards.sort()
    for i in range(len(cards) - 2):
        if is_consecutive(cards[i:i+3]):
            return 2  # Chi

    return 0  # Pass

# Initial cards
cards = list("PIZAUNOEWLWFQ")

# Rounds
rounds = [
    ('R', 'F'),
    ('R', 'I'),
    ('L', 'N'),
    ('S', 'L'),
    ('S', None)  # No card to remove in the final round
]

# Process each round
for new_card, remove_card in rounds:
    if remove_card:
        cards.remove(remove_card)
    result = determine_result(cards, new_card)
    cards.append(new_card)

# The result of the final round
print(result)