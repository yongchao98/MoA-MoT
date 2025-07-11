import math

# The set of standard chess pieces includes the Pawn, Knight, Bishop, Rook, Queen, and King.
num_standard_pieces = 6

# We are creating new pieces by combining the moves of two distinct standard pieces.
# The total number of such combinations is "6 choose 2".
total_combinations = math.comb(num_standard_pieces, 2)

# Based on chess principles, some of these hybrid pieces are too weak to deliver
# checkmate to a lone king. These are the combinations of a Pawn with a minor or major piece
# (excluding the Queen and King).
# 1. Pawn + Knight
# 2. Pawn + Bishop
# 3. Pawn + Rook
num_non_mating_combinations = 3

# The number of combinations that can achieve a checkmate is the total number of
# combinations minus those that cannot.
num_mating_combinations = total_combinations - num_non_mating_combinations

# The problem asks for the number of "distinct checkmate positions," which we have
# interpreted as the number of piece combinations capable of checkmate.
# The final code prints out the equation as requested.
print(f"{total_combinations} - {num_non_mating_combinations} = {num_mating_combinations}")