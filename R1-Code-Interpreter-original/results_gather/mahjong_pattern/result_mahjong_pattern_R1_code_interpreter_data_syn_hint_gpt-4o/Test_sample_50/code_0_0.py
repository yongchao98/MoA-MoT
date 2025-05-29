def check_peng_chi(cards, new_card):
    # Check for Peng
    if cards.count(new_card) >= 2:
        return 1  # Peng

    # Check for Chi
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return 2  # Chi

    return 0  # Pass

# Initial cards
cards = list("UMHYWCBRVQBJZ")

# Round 1
cards.remove('V')
result = check_peng_chi(cards, '[')

# Round 2
cards.remove('R')
result = check_peng_chi(cards, 'P')

# Round 3
cards.remove('Z')
result = check_peng_chi(cards, 'N')

# Round 4
cards.remove('Q')
result = check_peng_chi(cards, 'B')

# Round 5
result = check_peng_chi(cards, 'Y')

print(result)