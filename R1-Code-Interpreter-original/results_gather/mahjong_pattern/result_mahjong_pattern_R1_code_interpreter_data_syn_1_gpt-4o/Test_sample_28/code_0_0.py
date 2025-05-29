def check_peng(cards, new_card):
    # Check for "Peng"
    from collections import Counter
    card_count = Counter(cards)
    return card_count[new_card] >= 2

def check_chi(cards, new_card):
    # Check for "Chi"
    card_set = set(cards)
    new_card_ord = ord(new_card)
    # Check if there are two cards that can form a sequence with the new card
    return ((chr(new_card_ord - 1) in card_set and chr(new_card_ord - 2) in card_set) or
            (chr(new_card_ord - 1) in card_set and chr(new_card_ord + 1) in card_set) or
            (chr(new_card_ord + 1) in card_set and chr(new_card_ord + 2) in card_set))

def determine_result(cards, new_card):
    if check_peng(cards, new_card):
        return 1  # Peng
    elif check_chi(cards, new_card):
        return 2  # Chi
    else:
        return 0  # Pass

# Initial cards
cards = list("YZHDJJOTUUXLA")

# Round 1
cards.append('Z')
cards.remove('Y')
result1 = determine_result(cards, 'Z')

# Round 2
cards.append('Z')
cards.remove('Z')
result2 = determine_result(cards, 'Z')

# Round 3
cards.append('Z')
cards.remove('J')
result3 = determine_result(cards, 'Z')

# Round 4
cards.append('R')
cards.remove('A')
result4 = determine_result(cards, 'R')

# Round 5
cards.append('Z')
result5 = determine_result(cards, 'Z')

# Print the result of the final round
print(result5)