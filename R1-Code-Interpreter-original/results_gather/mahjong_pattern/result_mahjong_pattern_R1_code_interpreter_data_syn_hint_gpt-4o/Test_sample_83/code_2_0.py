def check_peng_chi(cards, new_card):
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
cards = list("OUFAEYRCKZYJW")

# Round 1
cards.append('E')
cards.remove('J')

# Round 2
cards.append('M')
cards.remove('Z')

# Round 3
cards.append('G')
cards.remove('C')

# Round 4
cards.append('E')
cards.remove('Y')

# Round 5
result_round_5 = check_peng_chi(cards, 'B')

print(result_round_5)