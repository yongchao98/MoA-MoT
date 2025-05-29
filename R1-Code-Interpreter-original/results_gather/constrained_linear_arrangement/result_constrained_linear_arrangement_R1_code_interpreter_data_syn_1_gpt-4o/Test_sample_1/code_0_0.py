# Deduce the opponent's cards based on the given rounds and results
opponent_cards = []

# Round 1: You played B (Wood), and you lost. Opponent must have played Fire (B).
opponent_cards.append('B')

# Round 2: You played C (Earth), and it was a draw. Opponent must have played Earth (C).
opponent_cards.append('C')

# Round 3: You played A (Wood), and it was a draw. Opponent must have played Wood (A).
opponent_cards.append('A')

# Print the deduced opponent's cards
print(opponent_cards)