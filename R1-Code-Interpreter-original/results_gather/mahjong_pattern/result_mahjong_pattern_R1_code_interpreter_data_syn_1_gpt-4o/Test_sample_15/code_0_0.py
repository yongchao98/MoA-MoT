def check_peng(cards, new_card):
    # Check for "Peng"
    card_count = {}
    for card in cards:
        if card in card_count:
            card_count[card] += 1
        else:
            card_count[card] = 1
    if card_count.get(new_card, 0) >= 2:
        return True
    return False

def check_chi(cards, new_card):
    # Check for "Chi"
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

# Final result for Round 4
print(result4)