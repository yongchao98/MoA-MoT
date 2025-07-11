# This script will print the best sequence of moves for the given chess puzzle.
# The puzzle describes a forced checkmate for Black.
# The best sequence is a classic smothered mate pattern.

# The moves are:
# 1... Qg1+ : Black sacrifices the Queen, forcing the Rook to capture.
# 2. Rxg1 : White's only legal move.
#    Nf2# : Black delivers checkmate with the Knight.

# The final sequence of moves is printed below.
print("The best sequence of moves is:")
print("1... Qg1+ 2. Rxg1 Nf2#")