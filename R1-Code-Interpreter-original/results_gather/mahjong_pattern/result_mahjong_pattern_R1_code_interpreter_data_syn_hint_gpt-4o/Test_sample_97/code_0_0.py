def is_consecutive(seq):
    """Check if the given sequence of letters is consecutive."""
    return ord(seq[1]) == ord(seq[0]) + 1 and ord(seq[2]) == ord(seq[1]) + 1

def determine_result(cards, new_card):
    """Determine the result based on the current cards and the new card."""
    # Check for Peng
    if cards.count(new_card) >= 2:
        return 1  # Peng

    # Check for Chi
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if is_consecutive(sorted_cards[i:i+3]):
            return 2  # Chi

    return 0  # Pass

# Initial cards
cards = list("IVAJUKTRNLQOH")

# Rounds
rounds = [
    ('S', 'N'),
    ('X', 'U'),
    ('B', 'Q'),
    ('C', 'S'),
    ('V', None)  # No card to remove in the final round
]

# Process each round
for new_card, remove_card in rounds:
    cards.append(new_card)
    if remove_card:
        cards.remove(remove_card)
    result = determine_result(cards, new_card)

# The result of the final round
print(result)