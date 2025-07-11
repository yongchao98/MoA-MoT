def solve_three_check_chess():
    """
    This function explains the solution to the three-check chess puzzle
    and prints the minimal number of moves for White to win.
    """

    explanation = """
The problem asks for the minimal number of moves for White to win, assuming optimal play from both sides. This requires finding a forced winning sequence for White while considering that Black will always choose the path that prolongs the game.

The winning line for White begins with the move 1. Bxd7+.

1.  White plays `Bxd7+`. This is the first check. (Checks: W=1, B=0)
    Black has two main replies: `Kxd7` or `Nxd7`.

    a) If Black plays `1... Kxd7`, White continues `2. O-O-O+` (Check 2). Black's king is forced to move, and on White's 3rd turn, White can deliver the third check and win (e.g., `2... Kc6 3. Qc4+` (Check 3, Win)). This line leads to a win for White in 3 moves.

    b) If Black plays `1... Nxd7`, White continues `2. Qb7+` (Check 2). Black must respond `2... Rb8`. White then plays the key move `3. Bxf6!`. Black cannot defend against all threats.
       - If Black plays `3... gxf6`, White delivers `4. Qxc8+` (Check 3, Win).
       - If Black plays `3... Rxb7`, White delivers `4. Bxe7+` (Check 3, Win).
       This line leads to a win for White in 4 moves.

Since Black plays optimally, they will choose option (b) to survive as long as possible, forcing the game to last 4 moves for White.

Therefore, the minimal number of moves for White to secure a win is 4.
    """
    print(explanation)

    # The final answer is an integer.
    # The "equation" might refer to the moves in the winning line.
    # Move 1: Bxd7+
    # Move 2: Qb7+
    # Move 3: Bxf6
    # Move 4: Qxc8+ (or Bxe7+)
    
    final_answer = 4
    
    print("Final Answer:")
    # As per the instruction "output each number in the final equation",
    # we will print the numbers of the moves leading to the win.
    print(1)
    print(2)
    print(3)
    print(4)

solve_three_check_chess()