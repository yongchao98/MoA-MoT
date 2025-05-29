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
cards = list("SNXQWCHBFYQZB")

# Round 1
cards.remove('W')
result1 = check_peng_chi(cards, '[')
cards.append('[')

# Round 2
cards.remove('Z')
result2 = check_peng_chi(cards, 'F')
cards.append('F')

# Round 3
cards.remove('B')
result3 = check_peng_chi(cards, 'U')
cards.append('U')

# Round 4
cards.remove('F')
result4 = check_peng_chi(cards, 'F')
cards.append('F')

# Round 5
result5 = check_peng_chi(cards, 'J')

print(result5)