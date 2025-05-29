# Initial cards
cards = list("MNZXZYOKWCFHU")

# Function to check for "Peng"
def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

# Function to check for "Chi"
def check_chi(cards, new_card):
    # Sort the cards and check for consecutive sequence
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return True
    return False

# Round 1
cards.append('[')
cards.remove('Z')
if check_peng(cards, '['):
    result_round_1 = 1
elif check_chi(cards, '['):
    result_round_1 = 2
else:
    result_round_1 = 0

# Round 2
cards.append('G')
cards.remove('Y')
if check_peng(cards, 'G'):
    result_round_2 = 1
elif check_chi(cards, 'G'):
    result_round_2 = 2
else:
    result_round_2 = 0

# Round 3
cards.append('I')
if check_peng(cards, 'I'):
    result_round_3 = 1
elif check_chi(cards, 'I'):
    result_round_3 = 2
else:
    result_round_3 = 0

# Output the result of the final round
print(result_round_3)