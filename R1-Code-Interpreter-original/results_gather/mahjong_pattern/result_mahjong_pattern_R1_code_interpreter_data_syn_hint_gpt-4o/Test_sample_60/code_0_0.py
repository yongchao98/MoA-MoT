def check_peng(cards, new_card):
    # Check if there are two identical cards and the new card is the same
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check if the new card can form a consecutive sequence with two other cards
    card_set = set(cards)
    new_card_ord = ord(new_card)
    return (chr(new_card_ord - 1) in card_set and chr(new_card_ord + 1) in card_set) or \
           (chr(new_card_ord - 2) in card_set and chr(new_card_ord - 1) in card_set) or \
           (chr(new_card_ord + 1) in card_set and chr(new_card_ord + 2) in card_set)

def determine_result(cards, new_card):
    if check_peng(cards, new_card):
        return 1  # Peng
    elif check_chi(cards, new_card):
        return 2  # Chi
    else:
        return 0  # Pass

# Initial cards
cards = list("GAWRMZMPOWJHV")

# Round 1
cards.append('J')
cards.remove('R')
determine_result(cards, 'J')

# Round 2
cards.append('W')
cards.remove('Z')
determine_result(cards, 'W')

# Round 3
cards.append('P')
cards.remove('P')
determine_result(cards, 'P')

# Round 4
cards.append('M')
cards.remove('O')
determine_result(cards, 'M')

# Round 5
cards.append('W')
result = determine_result(cards, 'W')

print(result)