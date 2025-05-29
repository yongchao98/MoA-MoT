def check_peng_or_chi(cards, new_card):
    # Check for Peng
    if cards.count(new_card) >= 2:
        return 1  # Peng

    # Check for Chi
    sorted_cards = sorted(cards + [new_card])
    for i in range(len(sorted_cards) - 2):
        if ord(sorted_cards[i+1]) == ord(sorted_cards[i]) + 1 and ord(sorted_cards[i+2]) == ord(sorted_cards[i]) + 2:
            return 2  # Chi

    return 0  # Pass

# Initial cards
cards = list("EJXXLNRZHAGDQ")

# Round 1
cards.append('I')
cards.remove('H')
result_round_1 = check_peng_or_chi(cards, 'I')

# Round 2
cards.append('F')
cards.remove('A')
result_round_2 = check_peng_or_chi(cards, 'F')

# Round 3
cards.append('K')
cards.remove('E')
result_round_3 = check_peng_or_chi(cards, 'K')

# Round 4
cards.append('G')
cards.remove('J')
result_round_4 = check_peng_or_chi(cards, 'G')

# Round 5
cards.append('R')
result_round_5 = check_peng_or_chi(cards, 'R')

print(result_round_5)