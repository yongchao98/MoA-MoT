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
cards = list("EJXXLNRZHAGDQ")

# Round 1
cards.remove('H')
cards.append('I')
result1 = determine_result(cards, 'I')

# Round 2
cards.remove('A')
cards.append('F')
result2 = determine_result(cards, 'F')

# Round 3
cards.remove('E')
cards.append('K')
result3 = determine_result(cards, 'K')

# Round 4
cards.remove('J')
cards.append('G')
result4 = determine_result(cards, 'G')

# Round 5
cards.append('R')
result5 = determine_result(cards, 'R')

print(result5)