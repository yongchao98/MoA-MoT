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
cards = list("TRSQOSEPUIHXD")

# Round 1
cards.remove('T')
result1 = check_peng_chi(cards, 'L')
cards.append('L')

# Round 2
cards.remove('R')
result2 = check_peng_chi(cards, 'G')
cards.append('G')

# Round 3
cards.remove('D')
result3 = check_peng_chi(cards, 'T')
cards.append('T')

# Round 4
cards.remove('P')
result4 = check_peng_chi(cards, 'J')
cards.append('J')

# Round 5
result5 = check_peng_chi(cards, 'U')

print(result5)