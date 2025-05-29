def check_peng(cards, new_card):
    # Check if there are two identical cards and the new card is the same
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check if the new card and two cards in hand can form a consecutive sequence
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
cards = list("WETWUCCPFNAIZ")

# Round 1
cards.append('W')
cards.remove('I')
result1 = determine_result(cards, 'W')

# Round 2
cards.append('N')
cards.remove('W')
result2 = determine_result(cards, 'N')

# Round 3
cards.append('C')
cards.remove('E')
result3 = determine_result(cards, 'C')

# Round 4
cards.append('G')
cards.remove('Z')
result4 = determine_result(cards, 'G')

# Round 5
cards.append('M')
result5 = determine_result(cards, 'M')

print(result5)