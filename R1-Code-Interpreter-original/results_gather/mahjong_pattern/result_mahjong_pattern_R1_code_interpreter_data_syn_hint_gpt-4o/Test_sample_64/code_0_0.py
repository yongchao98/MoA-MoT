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
cards = list("GJTQXOSWZHMPO")

# Round 1
cards.append('U')
cards.remove('M')
result_round_1 = check_peng_chi(cards, 'U')

# Round 2
cards.append('V')
cards.remove('G')
result_round_2 = check_peng_chi(cards, 'V')

# Round 3
cards.append('M')
cards.remove('U')
result_round_3 = check_peng_chi(cards, 'M')

# Round 4
cards.append('U')
cards.remove('W')
result_round_4 = check_peng_chi(cards, 'U')

# Round 5
cards.append('O')
result_round_5 = check_peng_chi(cards, 'O')

# Print the result of the final round
print(result_round_5)