def check_peng(cards, new_card):
    # Check for "Peng"
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check for "Chi"
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return True
    return False

def determine_result(cards, new_card):
    if check_peng(cards, new_card):
        return 1  # Peng
    elif check_chi(cards, new_card):
        return 2  # Chi
    else:
        return 0  # Pass

# Initial cards
cards = list("JAHXUVZSQWDCL")

# Round 1
cards.append('N')
cards.remove('S')
result1 = determine_result(cards, 'N')

# Round 2
cards.append('M')
cards.remove('J')
result2 = determine_result(cards, 'M')

# Round 3
cards.append('K')
cards.remove('C')
result3 = determine_result(cards, 'K')

# Round 4
cards.append('C')
cards.remove('L')
result4 = determine_result(cards, 'C')

# Round 5
cards.append('O')
result5 = determine_result(cards, 'O')

# Final result
print(result5)