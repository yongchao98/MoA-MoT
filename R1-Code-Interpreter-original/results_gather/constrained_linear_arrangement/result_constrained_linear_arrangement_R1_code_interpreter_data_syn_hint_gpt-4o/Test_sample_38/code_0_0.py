# Possible pieces for each round based on the deductions
round_1_possibilities = ['C', 'D', 'E', 'F', 'G', 'H']  # Smaller than B
round_2_possibilities = ['B', 'C', 'D', 'E', 'F', 'G', 'H']  # Smaller than A
round_3_piece = 'A'  # Opponent played A in Round 3

# We need to find a combination of pieces for Round 1 and Round 2
# that satisfies the constraints and does not repeat any piece.

# Iterate over possible combinations
for r1_piece in round_1_possibilities:
    for r2_piece in round_2_possibilities:
        # Ensure no repetition of pieces
        if r1_piece != r2_piece and r1_piece != round_3_piece and r2_piece != round_3_piece:
            print(f"Round 1: {r1_piece}, Round 2: {r2_piece}, Round 3: {round_3_piece}")