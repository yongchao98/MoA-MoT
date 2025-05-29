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
cards = list("NCBZWVTEMBIHU")

# Round 1
cards.append('J')
cards.remove('I')
result1 = determine_result(cards, 'J')

# Round 2
cards.append('X')
cards.remove('T')
result2 = determine_result(cards, 'X')

# Round 3
cards.append('D')
cards.remove('B')
result3 = determine_result(cards, 'D')

# Round 4
cards.append('H')
cards.remove('N')
result4 = determine_result(cards, 'H')

# Round 5
cards.append('H')
result5 = determine_result(cards, 'H')

print(result5)