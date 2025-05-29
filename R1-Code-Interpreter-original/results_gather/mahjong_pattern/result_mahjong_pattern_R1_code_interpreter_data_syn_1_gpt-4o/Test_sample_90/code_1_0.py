def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    sorted_cards = sorted(cards + [new_card])  # Convert new_card to a list
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
cards = list("NVSEXAXUITMHR")

# Round 1
cards.remove('M')
result1 = determine_result(cards, 'L')
cards.append('L')

# Round 2
cards.remove('S')
result2 = determine_result(cards, 'N')
cards.append('N')

# Round 3
cards.remove('H')
result3 = determine_result(cards, 'O')
cards.append('O')

# Round 4
cards.remove('O')
result4 = determine_result(cards, 'X')
cards.append('X')

# Round 5
result5 = determine_result(cards, 'V')

# Final result
print(result5)