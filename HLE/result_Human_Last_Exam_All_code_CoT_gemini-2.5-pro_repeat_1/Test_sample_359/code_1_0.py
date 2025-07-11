# The chess position is given by the FEN: 8/2k5/5pn1/1Pp1pNpp/3PP3/4K1B1/8/8 w - - 0 43
# The solution involves a forced sequence of moves for White to get the king to a central square.

# 1. White plays c6. This forces Black to respond.
print("1. c6")

# Black must capture the pawn with bxc6, otherwise White's king marches to d5 unopposed.
print("   ...bxc6")

# 2. White recaptures with dxc6. The path for the king is now clear.
print("2. dxc6")

# Black's best defense is to attack the c6-pawn with its king.
print("   ...Kb6")

# 3. White ignores the threat and continues the king march with Kd3.
print("3. Kd3")

# Black must capture the pawn.
print("   ...Kxc6")

# 4. White's king moves to the critical square Kc4.
print("4. Kc4")

# Black's king must move, for example to Kb6.
print("   ...Kb6")

# 5. White's king moves to d5, a central square, and wins the game.
print("5. Kd5")

# The number of moves White makes to win is 5.
winning_moves_count = 5
print(f"\nWhite can win in {winning_moves_count} moves.")
