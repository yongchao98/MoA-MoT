# Define piece values
queen_value = 9
rook_value = 5
pawn_value = 1

# Number of pieces in the constructed position
# 1. Mechanism Pieces
num_piston_rooks = 1
num_piston_queens = 5
num_caging_queens = 4

# 2. Pawn Wall
num_wall_pawns = 6

# 3. Filler Pieces
# Total squares = 64. 2 for kings, 2 for empty squares.
total_white_pieces = 60
num_filler_queens = total_white_pieces - (num_piston_rooks + num_piston_queens + num_caging_queens + num_wall_pawns)

# Calculate the total score
piston_score = (num_piston_rooks * rook_value) + (num_piston_queens * queen_value)
caging_score = num_caging_queens * queen_value
wall_score = num_wall_pawns * pawn_value
filler_score = num_filler_queens * queen_value

total_score = piston_score + caging_score + wall_score + filler_score

# Output the equation
total_queens = num_piston_queens + num_caging_queens + num_filler_queens

print(f"The calculation for the total material points is based on a board with {total_white_pieces} white pieces:")
print(f"({num_piston_rooks} Rook * {rook_value} pts) + ({total_queens} Queens * {queen_value} pts) + ({num_wall_pawns} Pawns * {pawn_value} pts)")
print(f"= ({num_piston_rooks * rook_value}) + ({total_queens * queen_value}) + ({num_wall_pawns * pawn_value})")
print(f"= {num_piston_rooks * rook_value + total_queens * queen_value + num_wall_pawns * pawn_value}")
print("\nFinal Answer:")
print(f"The greatest number of points is {total_score}.")
