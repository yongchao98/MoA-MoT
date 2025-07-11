# It's White's turn to move, and the goal is to find a checkmate in 2 moves.

# 1. White's first move: The Rook on f3 captures the pawn on f4, delivering a check.
#    1. Rf4+
#    This is a forcing move. The Black King on e4 is in check.
#    - The King cannot move to d3 (attacked by Bishop on g1).
#    - The King cannot move to f5 (attacked by Rook on h5).
#    - The check cannot be blocked.
#    - Black's only legal move is to capture the Rook.
#    1... Kxf4

# 2. White's second move: The Knight on e7 moves to d6, delivering checkmate.
#    2. Nfd6#
#    The Black King on f4 is now in check by the Knight on d6. This is checkmate because:
#    - The checking piece (Knight) cannot be captured.
#    - The check cannot be blocked.
#    - The King has no escape squares:
#      - e5 and f3 are attacked by the Knight on d6.
#      - d4 and d5 are attacked by the Queen on a4.
#      - e3 and d3 are attacked by the Bishop on g1.
#      - f5 and g5 are attacked by the Rook on h5.
#      - g4, g3, and g5 are attacked by the King on h4.
#      - g3 is also attacked by the Bishop on g2.

print("The 2-move checkmate is:")
print("1. Rf4+ Kxf4")
print("2. Nfd6#")