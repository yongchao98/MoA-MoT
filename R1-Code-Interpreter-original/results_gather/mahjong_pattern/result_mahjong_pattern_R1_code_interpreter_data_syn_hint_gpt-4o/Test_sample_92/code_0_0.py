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
cards = list("GRCPONYJVADUA")

# Round 1
cards.append('A')
cards.remove('C')
result_round_1 = check_peng_or_chi(cards, 'A')

# Round 2
cards.append('N')
cards.remove('Y')
result_round_2 = check_peng_or_chi(cards, 'N')

# Round 3
cards.append('W')
cards.remove('A')
result_round_3 = check_peng_or_chi(cards, 'W')

# Round 4
cards.append('Q')
cards.remove('D')
result_round_4 = check_peng_or_chi(cards, 'Q')

# Round 5
cards.append('S')
result_round_5 = check_peng_or_chi(cards, 'S')

# Output the result of the final round
print(result_round_5)