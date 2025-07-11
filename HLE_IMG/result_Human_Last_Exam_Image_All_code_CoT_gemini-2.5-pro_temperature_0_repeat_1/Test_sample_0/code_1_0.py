# This script prints the chess notation for the mate-in-2 sequence.
# The puzzle is for Black to move and mate in two, without moving the queens.

# The sequence is:
# 1. Black moves Knight to h2, capturing the pawn and delivering a check (Nxh2+).
#    White's only legal response is to capture the knight with the King (Kxh2).
# 2. Black moves Rook to f2, delivering checkmate (Rf2#).

print("The sequence for black is:")
print("1...Nxh2+")
print("2...Rf2#")