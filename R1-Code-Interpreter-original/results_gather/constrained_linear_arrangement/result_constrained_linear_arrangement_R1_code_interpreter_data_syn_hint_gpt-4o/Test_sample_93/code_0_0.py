# Possible pieces for each round based on the analysis
round_1_possibilities = {'C', 'D', 'F', 'G', 'H'}
round_3_possibilities = {'B', 'C', 'D', 'F', 'G'}
round_4_possibilities = {'A', 'B', 'C'}
round_5_possibilities = {'A', 'B'}

# E is already used in round 2
used_pieces = {'E'}

# Deduce the opponent's pieces
opponent_pieces = []

# Round 5 must be A or B
for piece in round_5_possibilities:
    if piece not in used_pieces:
        opponent_pieces.append(piece)
        used_pieces.add(piece)
        break

# Round 4 must be A, B, or C
for piece in round_4_possibilities:
    if piece not in used_pieces:
        opponent_pieces.append(piece)
        used_pieces.add(piece)
        break

# Round 3 must be B, C, D, F, or G
for piece in round_3_possibilities:
    if piece not in used_pieces:
        opponent_pieces.append(piece)
        used_pieces.add(piece)
        break

# Round 1 must be C, D, F, G, or H
for piece in round_1_possibilities:
    if piece not in used_pieces:
        opponent_pieces.append(piece)
        used_pieces.add(piece)
        break

# Round 2 is E, already known
opponent_pieces.append('E')

# Print the opponent's pieces in order of rounds
print(opponent_pieces)