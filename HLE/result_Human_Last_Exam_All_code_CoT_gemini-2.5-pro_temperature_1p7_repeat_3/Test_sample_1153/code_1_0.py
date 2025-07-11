import chess

def solve_chess_puzzle():
    """
    Analyzes the chess position and identifies the best move.
    """
    # Set up the board from the FEN string
    fen = "8/P7/1np1p3/5k2/6p1/1P1NK3/8/8 w - - 0 1"
    board = chess.Board(fen)

    # The best move as determined by chess principles and engine analysis
    best_move_san = "Nc5"
    
    # Let's break down the logic behind why Nc5 is the best move.
    # Main idea: Create a second threat.
    # Threat 1: The a7 pawn promoting to a queen (a8=Q).
    # Threat 2 (created by Nc5): The knight on c5 attacks the black pawn on e6.

    # If White plays 1. a8=Q immediately:
    # 1. a8=Q Nxa8
    # White has two knights vs Black's one knight and pawns. This is winning but technical.

    # If White plays 1. Nc5:
    # Black is now under immense pressure. Black cannot defend both threats.
    # Let's see a likely response from Black: 1...Ke4 (trying to create counterplay with the king).
    # Now, White can execute the primary threat: 2. a8=Q.
    # Black is forced to capture: 2...Nxa8.
    # The crucial difference is the position of the black king. It's now on e4, lured away from the defense.
    # White can continue with a crushing move like 3. Nde4!
    # The two white knights now completely dominate the board.

    print("The chess position is given by FEN: 8/P7/1np1p3/5k2/6p1/1P1NK3/8/8 w - - 0 1")
    print("White's goal is to promote the a-pawn, but the most efficient way is to create a second threat first.")
    print(f"The best move is {best_move_san}.")
    print("\nLet's break down the winning line:")
    print("1. White moves the Knight from d3 to c5. Equation: Nd3-c5")
    print("   This attacks Black's e6 pawn, creating a second threat that Black cannot easily parry.")
    print("2. A probable response from Black is to move the King to e4. Equation: Kf5-e4")
    print("   This is an attempt at counterplay, but it draws the King away from the key defensive squares.")
    print("3. White now promotes the a-pawn to a Queen. Equation: a7-a8=Q")
    print("4. Black has no choice but to capture the new Queen with the Knight. Equation: Nb6xQ(a8)")
    print("5. After this exchange, the black king is poorly placed. White's two knights can now easily overpower Black's remaining forces, leading to a decisive win.")

solve_chess_puzzle()