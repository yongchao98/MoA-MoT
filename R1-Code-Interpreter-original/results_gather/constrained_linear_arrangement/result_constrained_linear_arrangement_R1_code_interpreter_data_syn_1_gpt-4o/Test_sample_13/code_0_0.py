# Deduce the opponent's cards based on the given rounds and results
opponent_cards = []

# Round 1: We played Earth (E) and lost
# Earth (E) loses to Fire (B)
opponent_cards.append('B')

# Round 2: We played Metal (D) and it was a draw
# Draw means the opponent also played Metal (D)
opponent_cards.append('D')

# Round 3: We played Wood (A) and it was a draw
# Draw means the opponent also played Wood (A)
opponent_cards.append('A')

# Print the deduced opponent's cards
print(opponent_cards)