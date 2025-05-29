def check_peng_or_chi(cards, new_card):
    # Check for "Peng"
    if cards.count(new_card) >= 2:
        return 1  # Peng

    # Check for "Chi"
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return 2  # Chi

    return 0  # Pass

# Initial cards
cards = list("SDCBJJBLPELSN")

# Round 1
cards.append('S')
cards.remove('S')
result_round_1 = check_peng_or_chi(cards, 'S')

# Round 2
cards.append('L')
cards.remove('C')
result_round_2 = check_peng_or_chi(cards, 'L')

# Round 3
cards.append('S')
cards.remove('B')
result_round_3 = check_peng_or_chi(cards, 'S')

# Round 4
cards.append('H')
cards.remove('S')
result_round_4 = check_peng_or_chi(cards, 'H')

# Round 5
result_round_5 = check_peng_or_chi(cards, 'S')

print(result_round_5)