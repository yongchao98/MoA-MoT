# The chess position is r1b1k2r/ppp2ppp/8/2bBN3/3nn3/8/PPP2qPP/RNBQR2K b Qkq - 0 1
# It is black's turn to move and force a checkmate.
# The best sequence is a classic "Smothered Mate".

# Step 1: Black sacrifices the Queen to force the White Rook onto the g1 square.
black_move_1 = "1... Qg1+"

# Step 2: White is forced to capture the Queen, as it's the only legal move.
white_move_1 = "2. Rxg1"

# Step 3: Black's Knight delivers the checkmate. The White King is "smothered"
# by its own pieces and cannot escape the Knight's attack.
black_move_2 = "2... Nf2#"

print("The best sequence of moves is:")
print(f"Black moves: {black_move_1}")
print(f"White moves: {white_move_1}")
print(f"Black moves: {black_move_2}")