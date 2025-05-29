def check_peng(cards, new_card):
    # Check if there are two identical cards and the new card is the same
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check if the new card can form a consecutive sequence with any two cards
    cards_with_new = cards + [new_card]
    cards_with_new = sorted(set(cards_with_new))  # Sort and remove duplicates
    for i in range(len(cards_with_new) - 2):
        if ord(cards_with_new[i+1]) == ord(cards_with_new[i]) + 1 and ord(cards_with_new[i+2]) == ord(cards_with_new[i]) + 2:
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
cards = list("HZDFCXABNORCQ")

# Round 1
cards.remove('H')
cards.append('C')
result_round_1 = determine_result(cards, 'C')

# Round 2
cards.remove('C')
cards.append('O')
result_round_2 = determine_result(cards, 'O')

# Round 3
cards.remove('Q')
cards.append('F')
result_round_3 = determine_result(cards, 'F')

# Round 4
result_round_4 = determine_result(cards, 'C')

print(result_round_4)