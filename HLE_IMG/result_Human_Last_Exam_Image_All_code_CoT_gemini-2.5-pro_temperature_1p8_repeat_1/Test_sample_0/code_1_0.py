# This script prints the chess notation for the mate-in-2 sequence.
# 1. Black moves the Knight from f1 to g3, delivering a check (Ng3+).
#    White's only legal response is to move the King to h2.
# 2. Black moves the Rook from f8 to h8, delivering checkmate (Rh8#).
# The sequence is provided without the white move, as requested.

move1 = "Ng3+"
move2 = "Rh8#"
print(f"{move1}, {move2}")