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
cards = list("ENLVBRMAOZXKJ")

# Round 1
cards.append('P')
cards.remove('M')
result = check_peng_chi(cards, 'P')

# Round 2
cards.append('C')
cards.remove('L')
result = check_peng_chi(cards, 'C')

# Round 3
cards.append('L')
cards.remove('E')
result = check_peng_chi(cards, 'L')

# Round 4
cards.append('Q')
cards.remove('J')
result = check_peng_chi(cards, 'Q')

# Round 5
cards.append('R')
result = check_peng_chi(cards, 'R')

print(result)