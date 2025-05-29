def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
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
cards = list("DUACQSLJZSAWJ")

# Round 1
cards.append('J')
cards.remove('Q')
result1 = determine_result(cards, 'J')

# Round 2
cards.append('S')
cards.remove('A')
result2 = determine_result(cards, 'S')

# Round 3
cards.append('V')
cards.remove('J')
result3 = determine_result(cards, 'V')

# Round 4
cards.append('T')
cards.remove('Z')
result4 = determine_result(cards, 'T')

# Round 5
result5 = determine_result(cards, 'V')

print(result5)