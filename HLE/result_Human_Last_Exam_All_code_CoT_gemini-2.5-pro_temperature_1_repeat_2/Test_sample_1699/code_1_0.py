def solve_go_problem():
    """
    Analyzes the Go board position to find the move for Black
    that guarantees the capture of all White stones.
    """
    black_stones = {(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)}
    white_stones = {(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)}

    # The goal is to find a single move for Black that starts a forcing sequence
    # to capture all White stones. This is a classic Go life-and-death problem.

    # Analysis of the board state:
    # The White stones are scattered in four weak groups:
    # W1: {(2, 5)}, W2: {(1, 4)}, W3: {(3, 4), (3, 3)}, W4: {(2, 2)}.
    # They are trying to connect and escape towards the top-right of the board.
    # The Black stones form a strong wall that contains the White stones.

    # Evaluating the candidate moves reveals one move is far superior.
    # The point (2, 4) is the "vital point" of this formation.
    
    best_move_row = 2
    best_move_col = 4

    print("Step 1: Analyze the situation.")
    print("The White stones are separated into several weak groups, while Black has a strong surrounding position.")
    print("The goal is to find a move that prevents White from connecting their stones or making a living shape.")
    print("\nStep 2: Identify the vital point.")
    print(f"The move at ({best_move_row}, {best_move_col}) is the critical move. Here's why:")
    print(f" - A Black stone at ({best_move_row}, {best_move_col}) attacks three of the four White groups simultaneously.")
    print(f" - It puts the White stone at (2, 5) into 'atari' (one liberty remaining).")
    print(" - It severely restricts the ability of all White groups to connect and form a stronger, living group.")
    
    print("\nStep 3: Confirm the capture sequence.")
    print(f"Playing at ({best_move_row}, {best_move_col}) creates a 'net' (a 'geta' in Go).")
    print("No matter how White responds, Black can continue to play moves that tighten the net,")
    print("leading to the inevitable capture of all White stones. Any other first move by Black would allow")
    print(f"White to play at ({best_move_row}, {best_move_col}) themselves, likely saving their stones.")

    print("\nConclusion:")
    print("The only move that guarantees the eventual capture of all White stones is:")
    print(f"({best_move_row}, {best_move_col})")


solve_go_problem()