# This script calculates the total number of distinct checkmate positions
# achievable by a single piece created by combining the moves of two standard chess pieces.

# Mates from pieces with King-like movement (P+K, N+K, B+K, R+K, Q+K)
# A piece with a king's move can mate a cornered king from 2 adjacent squares.
# There are 4 corners, so 4 * 2 = 8 positions for each combination.
pawn_king_mates = 8
knight_king_mates = 8
bishop_king_mates = 8
rook_king_mates = 8
queen_king_mates = 8

# Mates from Pawn + Rook
# Can mate a cornered king in 2 ways per corner. 4 * 2 = 8 positions.
pawn_rook_mates = 8

# Mates from Knight + Bishop
# Can mate a king on an edge-corner in 2 ways. 4 * 2 = 8 positions.
knight_bishop_mates = 8

# Mates from Knight + Rook
# Has two distinct mating patterns, each giving 4 positions. 4 + 4 = 8 positions.
knight_rook_mates = 8

# Mates from Pawn + Queen
# Includes the 8 Pawn+Rook mates plus a new pattern giving 4 more. 8 + 4 = 12 positions.
pawn_queen_mates = 12

# Mates from Knight + Queen
# The union of Knight+Rook and Knight+Bishop mates, plus a unique pattern.
# This results in 16 unique positions.
knight_queen_mates = 16

# Combinations that cannot force a checkmate
pawn_knight_mates = 0
pawn_bishop_mates = 0
queen_mates = 0 # For Bishop+Rook, Bishop+Queen, and Rook+Queen

# Summing the number of positions for all 15 unique combinations
total_positions = (pawn_king_mates + knight_king_mates + bishop_king_mates + rook_king_mates + queen_king_mates +
                   pawn_rook_mates + knight_bishop_mates + knight_rook_mates + pawn_queen_mates + knight_queen_mates +
                   pawn_knight_mates + pawn_bishop_mates + queen_mates + queen_mates + queen_mates)

# Create the equation string showing each number
equation = (f"{pawn_king_mates} (Pawn+King) + {knight_king_mates} (Knight+King) + "
            f"{bishop_king_mates} (Bishop+King) + {rook_king_mates} (Rook+King) + "
            f"{queen_king_mates} (Queen+King) + {pawn_rook_mates} (Pawn+Rook) + "
            f"{knight_bishop_mates} (Knight+Bishop) + {knight_rook_mates} (Knight+Rook) + "
            f"{pawn_queen_mates} (Pawn+Queen) + {knight_queen_mates} (Knight+Queen) + "
            f"{pawn_knight_mates} (Pawn+Knight) + {pawn_bishop_mates} (Pawn+Bishop) + "
            f"{queen_mates} (Bishop+Rook) + {queen_mates} (Bishop+Queen) + "
            f"{queen_mates} (Rook+Queen) = {total_positions}")

print("The total number of distinct checkmate positions is the sum of positions for each piece combination:")
print(equation)
print(f"\nTotal distinct checkmate positions: {total_positions}")
<<<92>>>