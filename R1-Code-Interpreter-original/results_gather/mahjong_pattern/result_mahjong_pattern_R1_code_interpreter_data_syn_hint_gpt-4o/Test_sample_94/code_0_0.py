def check_peng(cards, new_card):
    # Check if there are two identical cards and the new card is the same
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check if two cards and the new card can form a consecutive sequence
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
cards = list("VFNFLXUIPRTCD")

# Round 1
cards.remove('N')
cards.append('E')
result1 = determine_result(cards, 'E')

# Round 2
cards.remove('U')
cards.append('G')
result2 = determine_result(cards, 'G')

# Round 3
cards.remove('P')
cards.append('C')
result3 = determine_result(cards, 'C')

# Round 4
cards.remove('G')
cards.append('F')
result4 = determine_result(cards, 'F')

# Round 5
cards.append('X')
result5 = determine_result(cards, 'X')

# Output the result of the final round
print(result5)