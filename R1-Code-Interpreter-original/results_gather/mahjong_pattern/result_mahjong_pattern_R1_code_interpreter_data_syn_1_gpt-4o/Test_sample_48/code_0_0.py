def check_peng(cards, new_card):
    # Check if there are two identical cards and the new card is the same
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check if the new card and any two cards can form a consecutive sequence
    cards_with_new = cards + [new_card]
    cards_with_new.sort()
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
cards = list("WBZGBEZNTIWFQ")

# Round 1
cards.remove('W')
cards.append('B')
result1 = determine_result(cards, 'B')

# Round 2
cards.remove('E')
cards.append('N')
result2 = determine_result(cards, 'N')

# Round 3
cards.remove('Z')
cards.append('U')
result3 = determine_result(cards, 'U')

# Round 4
cards.remove('T')
cards.append('X')
result4 = determine_result(cards, 'X')

# Round 5
result5 = determine_result(cards, 'Y')

print(result5)