def check_peng(cards, new_card):
    # Check if there are two identical cards and the new card is the same
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    # Check if the new card can form a consecutive sequence with any two cards
    cards.append(new_card)
    cards.sort()
    for i in range(len(cards) - 2):
        if ord(cards[i+1]) == ord(cards[i]) + 1 and ord(cards[i+2]) == ord(cards[i]) + 2:
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
cards = list("AWFDVZOZFPYCJ")

# Round 1
cards.append('U')
cards.remove('F')
result1 = determine_result(cards, 'U')

# Round 2
cards.append('K')
cards.remove('Y')
result2 = determine_result(cards, 'K')

# Round 3
cards.append('E')
cards.remove('Z')
result3 = determine_result(cards, 'E')

# Round 4
cards.append('E')
cards.remove('E')
result4 = determine_result(cards, 'E')

# Round 5
cards.append('G')
result5 = determine_result(cards, 'G')

print(result5)