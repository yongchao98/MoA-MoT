def check_peng(cards, new_card):
    # Check if there are two identical cards and the new card is the same
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check if the new card can form a consecutive sequence with any two cards
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
cards = list("BJMGVENDQLIYI")

# Round 1
cards.append('V')
cards.remove('I')
result_round_1 = determine_result(cards, 'V')

# Round 2
cards.append('N')
cards.remove('N')
result_round_2 = determine_result(cards, 'N')

# Round 3
cards.append('N')
cards.remove('M')
result_round_3 = determine_result(cards, 'N')

# Round 4
cards.append('N')
cards.remove('L')
result_round_4 = determine_result(cards, 'N')

# Round 5
cards.append('K')
result_round_5 = determine_result(cards, 'K')

# Print the result of the final round
print(result_round_5)