def check_peng(cards, new_card):
    # Check if there are two identical cards and the new card is the same
    card_count = {}
    for card in cards:
        if card in card_count:
            card_count[card] += 1
        else:
            card_count[card] = 1
    for card, count in card_count.items():
        if count >= 2 and card == new_card:
            return True
    return False

def check_chi(cards, new_card):
    # Check if two cards and the new card can form a consecutive sequence
    cards = sorted(cards + [new_card])
    for i in range(len(cards) - 2):
        if ord(cards[i+1]) == ord(cards[i]) + 1 and ord(cards[i+2]) == ord(cards[i]) + 2:
            return True
    return False

def determine_result(cards, new_card):
    if check_peng(cards, new_card):
        return 1  # Peng
    elif check_chi(cards, new_card):
        return 2  # Chi
    else:
        return 0  # Pass

# Initial cards
cards = list("VFNFLXUIPRTCD")

# Rounds
rounds = [
    ('E', 'N'),
    ('G', 'U'),
    ('C', 'P'),
    ('F', 'G'),
    ('X', None)  # In the last round, no card is removed
]

# Process each round
for new_card, remove_card in rounds:
    if remove_card:
        cards.remove(remove_card)
    result = determine_result(cards, new_card)
    cards.append(new_card)

# The result of the final round
print(result)