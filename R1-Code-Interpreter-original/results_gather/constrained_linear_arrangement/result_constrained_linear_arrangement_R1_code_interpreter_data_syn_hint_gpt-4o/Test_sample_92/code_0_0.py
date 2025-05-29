# Possible cards for each round based on the results
round_1_possibilities = {'B', 'C', 'D', 'E', 'F', 'G'}
round_2_possibilities = {'B', 'C', 'D', 'E', 'F'}
round_3_possibilities = {'B', 'C', 'D'}
round_4_possibilities = {'C', 'D', 'E', 'F', 'G', 'H'}
round_5_possibilities = {'E', 'F', 'G', 'H'}

# Find the intersection of possibilities for each round
round_1_card = round_1_possibilities.intersection(round_2_possibilities, round_3_possibilities)
round_2_card = round_2_possibilities.intersection(round_3_possibilities, round_4_possibilities)
round_3_card = round_3_possibilities.intersection(round_4_possibilities, round_5_possibilities)
round_4_card = round_4_possibilities.intersection(round_5_possibilities)
round_5_card = round_5_possibilities

# Determine the opponent's cards for each round
opponent_cards = [
    round_1_card.pop(),
    round_2_card.pop(),
    round_3_card.pop(),
    round_4_card.pop(),
    round_5_card.pop()
]

print(opponent_cards)