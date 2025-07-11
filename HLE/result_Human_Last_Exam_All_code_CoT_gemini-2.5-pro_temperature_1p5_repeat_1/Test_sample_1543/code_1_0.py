import sys

# This problem is a tactical chess puzzle that requires logical deduction rather than computation.
# The goal is to find the minimum number of moves for White to force a checkmate.

# The position from FEN 9k/5c1pb1/10/10/10/3Q6/PP5A2/K9 w - - 0 1:
# White: Ka1, Qd3, Pa2, Pb2, Ah2 (Archbishop: Bishop + Knight)
# Black: Kj8, Cf7 (Chancellor: Rook + Knight), Ph7, Bi7

# Let's analyze the shortest path to victory for White.

# White's most forcing move is 1. Qd8+, checking the Black king.
# Black has two primary ways to block this check. We must find a mate for both cases.

# Case 1: Black plays 1... Bf8
# 1. Qd8+   Bf8
# 2. Qxf8+  Ch8   (Forced, as Chancellor blocks the new check)
# 3. Qe7#   (This is checkmate. The Black Chancellor on h8 is pinned to the king on j8 by the White Archbishop on h2. The king's escape squares i7 and i8 are both covered by the White Queen on e7).
# This is a mate in 3 moves.

# Case 2: Black plays 1... Ch8
# This is Black's other optimal defense. White must also be able to force mate in 3 here.
# 1. Qd8+   Ch8
# 2. Ai4    (A powerful quiet move. White's Archbishop moves from h2 to i4)
# Black now has several possible moves, but White can mate on the next turn against all of them.
#   - If Black plays 2... Ki8 (moving the king to unpin the Chancellor), White plays 3. Qj7#, which is checkmate.
#   - If Black plays 2... h6 (pawn move), White plays 3. Qj6#, which is checkmate.
#   - If Black plays 2... Bh6 (bishop move, freeing i7), White plays 3. Qd7#, which is checkmate.

# Since White can force a checkmate in 3 moves regardless of Black's optimal defense,
# and a mate in 2 is not possible, the minimal number of moves to win is 3.

minimal_moves_to_win = 3

print(minimal_moves_to_win)
