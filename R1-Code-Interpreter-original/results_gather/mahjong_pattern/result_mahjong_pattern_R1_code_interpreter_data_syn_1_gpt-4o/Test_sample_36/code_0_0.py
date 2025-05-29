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
cards = list("SRNVLXKWYKRPB")

# Round 1
cards.append('T')
cards.remove('N')
result1 = determine_result(cards, 'T')

# Round 2
cards.append('Z')
cards.remove('R')
result2 = determine_result(cards, 'Z')

# Round 3
cards.append('V')
cards.remove('W')
result3 = determine_result(cards, 'V')

# Round 4
cards.append('R')
cards.remove('T')
result4 = determine_result(cards, 'R')

# Round 5
result5 = determine_result(cards, 'R')

print(result5)