# It's a mate in 2 puzzle for Black.
# Black must deliver checkmate in two moves, without moving the queens.

# 1. Black's first move is Knight from d1 to f2, delivering a check.
#    - Notation: Nf2+
#    - This forces the White King from h1 to g1, as h2 is controlled by Black's Queen on b1.

# 2. Black's second move is Knight from f2 to h3, delivering checkmate.
#    - Notation: Nh3#
#    - The White King on g1 is in check.
#    - All escape squares are covered:
#      - f1 by the Rook on f8.
#      - f2 by the Knight on h3.
#      - g2 by the Queen on a2.
#      - h1 by the Queen on b1 and the Knight on h3.
#      - h2 by the Queen on b1.
#    - The Knight on h3 cannot be captured. This is checkmate.

move1 = "Nf2+"
move2 = "Nh3#"

print(f"{move1}, {move2}")