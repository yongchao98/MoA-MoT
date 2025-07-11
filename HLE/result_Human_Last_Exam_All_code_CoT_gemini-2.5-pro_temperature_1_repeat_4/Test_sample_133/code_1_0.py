def find_optimal_moves():
    """
    Analyzes the specific Connect 4 board state to find the optimal moves for 'O'.

    An optimal move is one that leads to a win as fast as possible.

    The analysis is based on the following key board features:
    - O-Threat: A diagonal of three 'O's (c6, d5, e4) threatens a win at f3.
               This requires the f4 spot to be filled first.
    - X-Threat: A diagonal of two 'X's (d6, e5) threatens to become a line of three
               if X plays at f4, which would then threaten a win for X at g3.

    We will analyze O's legal moves (c4, f4, etc.) to see which ones secure a win.
    """

    # Analysis of move 'c4':
    # 1. O plays c4. This creates a horizontal line of three O's (c4, d4, e4).
    # 2. X MUST block this threat by playing at f4. If X plays anywhere else, O wins by playing f4.
    # 3. After X plays f4, the spot f3 becomes available.
    # 4. O can now play f3 to complete the diagonal (c6, d5, e4, f3) and win.
    # This is a forced win for O on their next turn.
    move1 = "c4"

    # Analysis of move 'f4':
    # 1. O plays f4. This move simultaneously:
    #    a) Blocks X's potential threat on the d6-e5 diagonal.
    #    b) Creates a horizontal line of three O's (d4, e4, f4).
    # 2. X MUST block this horizontal threat by playing at either c4 or g4.
    # 3. Whichever move X makes, the spot f3 is now available for O.
    # 4. O can play f3 to complete the diagonal (c6, d5, e4, f3) and win.
    # This is also a forced win for O on their next turn.
    move2 = "f4"
    
    # Other moves are suboptimal because they either lead to a slower win or allow X
    # to play f4, creating a winning threat for X that leads to O losing.
    # Since both 'c4' and 'f4' lead to a win in the same number of turns (2),
    # they are both optimal.

    optimal_moves = sorted([move1, move2])
    
    print(", ".join(optimal_moves))

find_optimal_moves()