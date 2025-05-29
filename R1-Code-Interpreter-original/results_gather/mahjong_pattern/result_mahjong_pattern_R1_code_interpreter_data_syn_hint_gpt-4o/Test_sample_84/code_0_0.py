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
cards = list("BOEQDLYDPBWNP")

# Round 1
cards.append('Z')
cards.remove('P')
result1 = check_peng_chi(cards, 'Z')

# Round 2
cards.append('C')
cards.remove('O')
result2 = check_peng_chi(cards, 'C')

# Round 3
cards.append('P')
cards.remove('B')
result3 = check_peng_chi(cards, 'P')

# Round 4
cards.append('P')
cards.remove('W')
result4 = check_peng_chi(cards, 'P')

# Round 5
cards.append('P')
result5 = check_peng_chi(cards, 'P')

# Print the result of the final round
print(result5)