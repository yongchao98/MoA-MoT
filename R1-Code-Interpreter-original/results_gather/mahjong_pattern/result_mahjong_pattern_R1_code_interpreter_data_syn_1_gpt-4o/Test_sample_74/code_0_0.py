def check_peng_or_chi(cards, new_card):
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
cards = list("AVRNHDROSGQYP")

# Round 1
cards.remove('H')
cards.append('I')
result_round_1 = check_peng_or_chi(cards, 'I')

# Round 2
cards.remove('V')
cards.append('T')
result_round_2 = check_peng_or_chi(cards, 'T')

# Round 3
cards.remove('Y')
cards.append('R')
result_round_3 = check_peng_or_chi(cards, 'R')

# Round 4
cards.remove('P')
cards.append('S')
result_round_4 = check_peng_or_chi(cards, 'S')

# Round 5
cards.append('M')
result_round_5 = check_peng_or_chi(cards, 'M')

print(result_round_5)