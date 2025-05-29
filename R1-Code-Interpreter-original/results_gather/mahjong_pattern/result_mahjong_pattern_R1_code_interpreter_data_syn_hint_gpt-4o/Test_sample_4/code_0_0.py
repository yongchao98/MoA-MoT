# Initial set of cards
cards = list("ASCQPGYIKNVBD")

# Round 1: Add 'I', remove 'D'
cards.append('I')
cards.remove('D')

# Check for "Peng" after Round 1
def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

# Round 2: Add 'V'
new_card = 'V'
cards.append(new_card)

# Check for "Peng" after Round 2
result = check_peng(cards, new_card)

# Output the result for the final round
print(1 if result else 0)