def check_peng(cards, new_card):
    # Check if there are two identical cards in hand and the new card is the same
    from collections import Counter
    card_count = Counter(cards)
    for card, count in card_count.items():
        if count >= 2 and card == new_card:
            return True
    return False

def check_chi(cards, new_card):
    # Check if the new card can form a consecutive sequence with any two cards
    card_values = [ord(card) for card in cards]
    new_card_value = ord(new_card)
    card_values.append(new_card_value)
    card_values = sorted(set(card_values))  # Sort and remove duplicates
    for i in range(len(card_values) - 2):
        if card_values[i+1] == card_values[i] + 1 and card_values[i+2] == card_values[i] + 2:
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
cards = list("AIFPSFSLWOCQO")

# Round 1
cards.append('E')
cards.remove('S')
result1 = determine_result(cards, 'E')

# Round 2
cards.append('M')
cards.remove('E')
result2 = determine_result(cards, 'M')

# Round 3
cards.append('W')
cards.remove('O')
result3 = determine_result(cards, 'W')

# Round 4
cards.append('W')
cards.remove('Q')
result4 = determine_result(cards, 'W')

# Round 5
cards.append('W')
result5 = determine_result(cards, 'W')

print(result5)