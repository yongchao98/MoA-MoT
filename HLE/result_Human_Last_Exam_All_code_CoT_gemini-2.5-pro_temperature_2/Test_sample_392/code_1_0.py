# The chess position is r1b1k2r/ppp2ppp/8/2bBN3/3nn3/8/PPP2qPP/RNBQR2K b Qkq - 0 1.
# It is black's turn to move and force a checkmate.
# The optimal sequence is a classic smothered mate pattern.

# Step 1: Black sacrifices the Queen to force the White Rook onto the g1 square.
move_1_black_number = "1..."
move_1_black_action = "Qg1+"

# Step 2: White has no choice but to capture the Queen with the Rook.
move_2_white_number = "2."
move_2_white_action = "Rxg1"

# Step 3: With the White King trapped by its own pieces, Black's Knight delivers the checkmate.
move_2_black_action = "Nf2#"

print("The best sequence for Black to force checkmate is:")
print(move_1_black_number, move_1_black_action, move_2_white_number, move_2_white_action, move_2_black_action)