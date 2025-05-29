def check_peng(cards, new_card):
    # Check if there are two identical cards in hand and the new card is the same
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check if the new card can form a consecutive sequence with any two cards in hand
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
cards = list("SNXQWCHBFYQZB")

# Round 1
cards.remove('W')
cards.append('[')
result1 = determine_result(cards, '[')

# Round 2
cards.remove('Z')
cards.append('F')
result2 = determine_result(cards, 'F')

# Round 3
cards.remove('B')
cards.append('U')
result3 = determine_result(cards, 'U')

# Round 4
cards.remove('F')
cards.append('F')
result4 = determine_result(cards, 'F')

# Round 5
cards.append('J')
result5 = determine_result(cards, 'J')

print(result5)