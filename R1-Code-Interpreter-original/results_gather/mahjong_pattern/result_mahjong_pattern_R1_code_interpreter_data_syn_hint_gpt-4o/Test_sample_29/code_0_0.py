def check_peng(cards, new_card):
    # Check if there are two identical cards and the new card is the same
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check if the new card can form a consecutive sequence with two other cards
    card_set = set(cards)
    new_card_ord = ord(new_card)
    return (chr(new_card_ord - 1) in card_set and chr(new_card_ord - 2) in card_set) or \
           (chr(new_card_ord - 1) in card_set and chr(new_card_ord + 1) in card_set) or \
           (chr(new_card_ord + 1) in card_set and chr(new_card_ord + 2) in card_set)

def determine_result(cards, new_card):
    if check_peng(cards, new_card):
        return 1  # Peng
    elif check_chi(cards, new_card):
        return 2  # Chi
    else:
        return 0  # Pass

# Initial cards
cards = list("ALUTPJGKOSGPJ")

# Round 1
cards.append('W')
cards.remove('J')
result1 = determine_result(cards, 'W')

# Round 2
cards.append('M')
cards.remove('G')
result2 = determine_result(cards, 'M')

# Round 3
cards.append('M')
cards.remove('O')
result3 = determine_result(cards, 'M')

# Round 4
cards.append('B')
cards.remove('A')
result4 = determine_result(cards, 'B')

# Round 5
cards.append('L')
result5 = determine_result(cards, 'L')

# Print the result of the final round
print(result5)