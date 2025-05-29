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
cards = list("PICDHPZEJGDMW")

# Round 1
cards.append('F')
cards.remove('M')
result1 = determine_result(cards, 'F')

# Round 2
cards.append('E')
cards.remove('Z')
result2 = determine_result(cards, 'E')

# Round 3
cards.append('E')
cards.remove('C')
result3 = determine_result(cards, 'E')

# Round 4
cards.append('K')
result4 = determine_result(cards, 'K')

print(result4)