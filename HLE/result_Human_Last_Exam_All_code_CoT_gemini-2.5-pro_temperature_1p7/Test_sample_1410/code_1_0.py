import sys

def solve_chess_puzzle():
    """
    This function explains and prints the solution to the chess puzzle.
    """
    winning_line = [
        "1. Qxe1 Rxe1",
        "2. Bxf7+ Kh8",
        "3. Nf7+ Kg8",
        "4. Nh6++ Kh8",
        "5. Rf8#"
    ]

    print("White can win in 5 moves.")
    print("The winning sequence is as follows:")
    
    # This loop prints each move, satisfying the "output each number" instruction
    for move in winning_line:
        # For Python 2 compatibility, we need to ensure this works.
        # However, the environment is likely Python 3.
        # The following line will work in both Python 2 and 3 if this script were executed.
        sys.stdout.write(move + "\n")
        
    print("\nThus, the answer is 5.")

solve_chess_puzzle()
