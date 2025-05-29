def check_peng(cards, new_card):
    # Check for "Peng"
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check for "Chi"
    card_set = set(cards)
    new_card_ord = ord(new_card)
    # Check for sequences like ABC, BCD, etc.
    if (chr(new_card_ord - 1) in card_set and chr(new_card_ord - 2) in card_set) or \
       (chr(new_card_ord - 1) in card_set and chr(new_card_ord + 1) in card_set) or \
       (chr(new_card_ord + 1) in card_set and chr(new_card_ord + 2) in card_set):
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
cards = list("OQVJEDCGCGBPO")

# Round 1
cards.append('T')
cards.remove('B')
result1 = determine_result(cards, 'T')

# Round 2
cards.append('R')
cards.remove('E')
result2 = determine_result(cards, 'R')

# Round 3
cards.append('O')
cards.remove('R')
result3 = determine_result(cards, 'O')

# Round 4
cards.append('S')
cards.remove('C')
result4 = determine_result(cards, 'S')

# Round 5
cards.append('G')
result5 = determine_result(cards, 'G')

print(result5)