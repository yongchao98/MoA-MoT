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
cards = list("UXCIWSLGOYDBE")

# Round 1
cards.append('D')
cards.remove('Y')
result_round_1 = check_peng_or_chi(cards, 'D')

# Round 2
cards.append('W')
cards.remove('O')
result_round_2 = check_peng_or_chi(cards, 'W')

# Round 3
cards.append('W')
cards.remove('I')
result_round_3 = check_peng_or_chi(cards, 'W')

# Round 4
cards.append('Y')
cards.remove('S')
result_round_4 = check_peng_or_chi(cards, 'Y')

# Round 5
cards.append('W')
result_round_5 = check_peng_or_chi(cards, 'W')

print(result_round_5)