def is_consecutive(seq):
    return all(ord(seq[i]) + 1 == ord(seq[i + 1]) for i in range(len(seq) - 1))

def check_chi(cards, new_card):
    cards.append(new_card)
    cards.sort()
    for i in range(len(cards) - 2):
        if is_consecutive(cards[i:i+3]):
            return True
    return False

# Initial cards after all rounds
cards = list("SQOSEUIHXLTGJ")
new_card = 'U'

# Check for Chi
result = 2 if check_chi(cards, new_card) else 0
print(result)