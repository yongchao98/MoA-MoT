import sys
# This problem asks to find the forced checkmate sequence for Black.
# After analyzing all the provided options, the correct sequence is a famous pattern
# known as a "smothered mate" or "Philidor's Legacy".

# Here is the breakdown of the correct sequence (Option C):

# 1... Qg1+
# Black begins with a queen sacrifice by moving the queen from f2 to g1,
# delivering a check to the white king on h1. White cannot capture the queen with the
# king because the g1 square is also attacked by the black knight on e4.
# This forces white's rook on f1 to capture the queen.

# 2. Rxg1
# White's rook captures the black queen. This move is forced. Now, the white king
# is trapped on h1, blocked in by its own rook on g1 and its own pawn on g2.

# 2... Nf2#
# Black's knight on d4 moves to f2, delivering checkmate. The white king cannot
# move to any escape squares, the checking knight cannot be captured, and the
# check cannot be blocked. This concludes the game.

# The other options are incorrect due to illegal moves or sequences that do not result in a forced mate.

# Let's print the correct sequence in Algebraic Classic Notation.
move1_black = "Qg1+"
move2_white = "Rxg1"
move3_black = "Nf2#"

print(f"1... {move1_black} 2. {move2_white} {move3_black}")

# Flush the output to ensure it's displayed immediately.
sys.stdout.flush()