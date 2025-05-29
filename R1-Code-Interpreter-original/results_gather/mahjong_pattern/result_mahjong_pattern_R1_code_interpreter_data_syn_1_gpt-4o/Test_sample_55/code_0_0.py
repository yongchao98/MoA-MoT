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
cards = list("PARMYQEXCGSRN")

# Round 1
cards.remove('Y')
cards.append('I')
result_round_1 = determine_result(cards, 'I')

# Round 2
cards.remove('R')
cards.append('Q')
result_round_2 = determine_result(cards, 'Q')

# Round 3
cards.remove('R')
cards.append('R')
result_round_3 = determine_result(cards, 'R')

# Round 4
cards.remove('I')
cards.append('N')
result_round_4 = determine_result(cards, 'N')

# Round 5
cards.append('G')
result_round_5 = determine_result(cards, 'G')

# Final result for Round 5
print(result_round_5)