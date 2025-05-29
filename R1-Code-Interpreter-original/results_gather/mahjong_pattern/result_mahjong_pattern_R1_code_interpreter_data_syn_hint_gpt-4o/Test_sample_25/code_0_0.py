def check_peng(cards, new_card):
    return cards.count(new_card) >= 2

def check_chi(cards, new_card):
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return True
    return False

# Initial cards
cards = list("JTXDYVEGYPKZ")

# Round 5: Add H
new_card = 'H'

# Check for Peng
if check_peng(cards, new_card):
    result = 1
# Check for Chi
elif check_chi(cards, new_card):
    result = 2
else:
    result = 0

print(result)