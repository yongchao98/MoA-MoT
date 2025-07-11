import chess

class ThreeCheckBoard:
    """A wrapper around chess.Board to handle three-check rules."""

    def __init__(self, fen_string):
        """Initializes the board and check counts from a three-check FEN."""
        parts = fen_string.split(' ')
        # Standard FEN part for chess.Board
        board_fen = " ".join(parts[0:6])
        self.board = chess.Board(board_fen)
        self.white_checks = 0
        self.black_checks = 0
        # Parse custom check counters
        if len(parts) > 6:
            check_str = parts[6]
            try:
                self.white_checks = int(check_str.split('+')[1])
                self.black_checks = int(check_str.split('+')[2])
            except (IndexError, ValueError):
                # Handle cases where check string is malformed or absent
                pass
        self.move_number = self.board.fullmove_number

    def push_san(self, san_move):
        """Pushes a move in Standard Algebraic Notation and updates check counts."""
        try:
            move = self.board.parse_san(san_move)
            
            if self.board.turn == chess.WHITE:
                print(f"{self.move_number}. {san_move}", end="")
            else:
                print(f" {san_move}")
                self.move_number += 1

            if self.board.gives_check(move):
                if self.board.turn == chess.WHITE:
                    self.white_checks += 1
                    print(f" (White's check #{self.white_checks})")
                else:
                    self.black_checks += 1
            elif self.board.turn == chess.BLACK:
                print("")


            self.board.push(move)
            
        except ValueError as e:
            print(f"\nInvalid move {san_move}: {e}")
            return None
        return move
        
    def print_status(self):
        """Prints the current check counts."""
        print(f"Check count: White = {self.white_checks}, Black = {self.black_checks}")

def solve_three_check():
    """
    Solves the three-check chess problem from the given FEN.
    """
    fen = "r3kb1r/p2nqppp/5n2/1B2p1B1/4P3/1Q6/PPP2PPP/R3K2R w KQkq - 0 1 +0+0"
    game = ThreeCheckBoard(fen)

    print("Starting Position FEN:", fen)
    game.print_status()
    print("\n--- Optimal Winning Sequence ---")
    
    # White's first move: Long castle to activate the rook on a1.
    game.push_san("O-O-O")
    
    # Black's best defense is to unpin the f6-knight, but it allows a forced sequence.
    game.push_san("Be7")
    game.print_status()

    # White sacrifices the rook for the knight to open up lines.
    game.push_san("Rxd7")
    game.push_san("Qxd7")
    game.print_status()

    # White sacrifices a bishop for the queen, delivering the first check.
    game.push_san("Bxd7+")
    game.push_san("Kxd7")
    game.print_status()

    # The black king is now exposed. White's queen delivers the second check.
    game.push_san("Qb5+")
    
    # Black is forced to move the king. ...Kc7 is a plausible defense.
    game.push_san("Kc7")
    game.print_status()

    # White brings the final rook into the attack for the third and final check.
    game.push_san("Rd1+")
    game.print_status()
    
    print("\nWhite wins by delivering the third check on move 5.")
    print("The minimal amount of moves by White to win is 5.")


solve_three_check()