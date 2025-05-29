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
cards = list("JAHXUVZSQWDCL")

# Round 1
cards.append('N')
cards.remove('S')
result1 = check_peng_chi(cards, 'N')

# Round 2
cards.append('M')
cards.remove('J')
result2 = check_peng_chi(cards, 'M')

# Round 3
cards.append('K')
cards.remove('C')
result3 = check_peng_chi(cards, 'K')

# Round 4
cards.append('C')
cards.remove('L')
result4 = check_peng_chi(cards, 'C')

# Round 5
cards.append('O')
result5 = check_peng_chi(cards, 'O')

print(result5)