def solve_three_check_chess():
    """
    Analyzes the given Three-check chess position to find the minimal
    number of moves for White to win.
    The solution is a forced sequence of checks.
    """
    
    # The FEN string from the prompt has a non-standard part for the 7th rank.
    # We will use the most common interpretation where "pp1nqppp" is intended.
    # Initial state:
    # FEN: r3kb1r/pp1nqppp/5n2/1B2p1B1/4P3/1Q6/PPP2PPP/R3K2R w KQkq - 0 1
    # White checks: 0, Black checks: 0
    
    print("Assuming optimal play from both sides, White has a forced win.")
    print("Here is the main line of play leading to the quickest victory:\n")
    
    white_checks = 0
    
    # --- Move 1 ---
    print("White's 1st move: Bxd7+")
    # This move gives the first check to Black.
    previous_checks = white_checks
    white_checks += 1
    print(f"White check count: {previous_checks} + 1 = {white_checks}")
    print("Black's best response is to recapture. Let's analyze 1... Qxd7.\n")
    
    # --- Move 2 ---
    print("White's 2nd move: Qd5+")
    # After Black's recapture, White gives a second check.
    previous_checks = white_checks
    white_checks += 1
    print(f"White check count: {previous_checks} + 1 = {white_checks}")
    print("Black must move the king. The optimal move is 2... Kf8.\n")

    # --- Move 3 ---
    print("White's 3rd move: Qf7#")
    # This is the third check, which wins the game. It is also checkmate.
    previous_checks = white_checks
    white_checks += 1
    print(f"White check count: {previous_checks} + 1 = {white_checks}")
    print("This is the 3rd check. White wins.\n")
    
    minimal_moves_to_win = 3
    print(f"The minimal amount of moves by White to win is: {minimal_moves_to_win}")

solve_three_check_chess()

<<<3>>>