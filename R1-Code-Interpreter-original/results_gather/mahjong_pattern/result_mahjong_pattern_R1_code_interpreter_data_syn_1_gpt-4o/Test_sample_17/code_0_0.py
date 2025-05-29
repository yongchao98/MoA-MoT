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
cards = list("SCABRLDAVXJHZ")

# Round 1
cards.append('K')
cards.remove('H')
result_round_1 = check_peng_chi(cards, 'K')

# Round 2
cards.append('M')
cards.remove('L')
result_round_2 = check_peng_chi(cards, 'M')

# Round 3
cards.append('C')
cards.remove('R')
result_round_3 = check_peng_chi(cards, 'C')

# Round 4
cards.append('C')
result_round_4 = check_peng_chi(cards, 'C')

# Output the result of the final round
print(result_round_4)