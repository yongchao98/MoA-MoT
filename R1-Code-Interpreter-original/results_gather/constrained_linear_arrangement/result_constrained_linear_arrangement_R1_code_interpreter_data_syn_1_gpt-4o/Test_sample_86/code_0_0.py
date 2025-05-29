# Possible pieces for each round based on the analysis
rounds = {
    1: ['A', 'B', 'C'],  # Larger than D
    2: ['B', 'D', 'E', 'F', 'G', 'H'],  # Smaller than A
    3: ['C'],  # Draw with C
    4: ['G', 'H'],  # Smaller than F
    5: ['E', 'F', 'D', 'B', 'A']  # Larger than G
}

# Deduce the opponent's pieces
opponent_pieces = [''] * 5

# Round 3 is known
opponent_pieces[2] = 'C'

# Round 4: G or H
if 'G' in rounds[4]:
    opponent_pieces[3] = 'G'
else:
    opponent_pieces[3] = 'H'

# Round 1: A, B, or C (C is used in Round 3)
if 'A' in rounds[1]:
    opponent_pieces[0] = 'A'
elif 'B' in rounds[1]:
    opponent_pieces[0] = 'B'

# Round 2: B, D, E, F, G, or H (C is used in Round 3)
for piece in rounds[2]:
    if piece not in opponent_pieces:
        opponent_pieces[1] = piece
        break

# Round 5: E, F, D, B, or A (C is used in Round 3)
for piece in rounds[5]:
    if piece not in opponent_pieces:
        opponent_pieces[4] = piece
        break

print(opponent_pieces)