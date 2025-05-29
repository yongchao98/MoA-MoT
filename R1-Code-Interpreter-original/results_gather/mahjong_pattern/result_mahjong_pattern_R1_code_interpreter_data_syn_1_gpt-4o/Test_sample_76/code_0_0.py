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
cards = list("IKTNTHZGBSFRF")

# Round 1
cards.append('F')
cards.remove('S')
result_round_1 = determine_result(cards, 'F')

# Round 2
cards.append('F')
cards.remove('F')
result_round_2 = determine_result(cards, 'F')

# Round 3
cards.append('F')
cards.remove('F')
result_round_3 = determine_result(cards, 'F')

# Round 4
cards.append('T')
cards.remove('H')
result_round_4 = determine_result(cards, 'T')

# Round 5
cards.append('A')
result_round_5 = determine_result(cards, 'A')

print(result_round_5)