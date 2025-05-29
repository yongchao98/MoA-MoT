def check_peng_chi(cards, new_card):
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
cards = list("SAVCKNWPBXQER")

# Round 1
cards.append('Y')
cards.remove('P')
result_round_1 = check_peng_chi(cards, 'Y')

# Round 2
cards.append('X')
cards.remove('S')
result_round_2 = check_peng_chi(cards, 'X')

# Round 3
cards.append('A')
cards.remove('W')
result_round_3 = check_peng_chi(cards, 'A')

# Round 4
cards.append('C')
cards.remove('V')
result_round_4 = check_peng_chi(cards, 'C')

# Round 5
cards.append('I')
result_round_5 = check_peng_chi(cards, 'I')

# Output the result of the final round
print(result_round_5)