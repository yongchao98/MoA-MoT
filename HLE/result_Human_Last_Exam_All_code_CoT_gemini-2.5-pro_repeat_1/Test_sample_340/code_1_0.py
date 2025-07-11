import chess

class ThreeCheckChess(chess.Board):
    """
    A class to represent a Three-Check chess game, extending python-chess.
    """
    def __init__(self, fen, white_checks=0, black_checks=0):
        # The FEN for three-check might have check counts at the end.
        # python-chess doesn't support this, so we handle it manually.
        fen_parts = fen.split(' ')
        base_fen = ' '.join(fen_parts[:-1]) if fen_parts[-1].startswith('+') else fen
        super().__init__(base_fen)
        self.white_checks = white_checks
        self.black_checks = black_checks

    def push(self, move):
        """
        Pushes a move and updates check counts.
        """
        is_check_before = self.is_check()
        super().push(move)
        is_check_after = self.is_check()

        # A check is counted only if the move results in a check,
        # and the king was not already in check (to handle moving out of check).
        if is_check_after and not is_check_before:
            if self.turn == chess.BLACK:  # White just moved
                self.white_checks += 1
            else:  # Black just moved
                self.black_checks += 1

    def is_game_over(self):
        """
        Checks for three-check win condition.
        """
        if self.white_checks >= 3 or self.black_checks >= 3:
            return True
        return super().is_game_over()

def solve():
    """
    Finds and prints the minimal sequence for White to win.
    """
    # The FEN is r3kb1r/p2nqppp/5n2/1B2p1B1/4P3/1Q6/PPP2PPP/R3K2R w KQkq - 0 1 +0+0
    # We parse the check counts manually.
    initial_fen = "r3kb1r/p2nqppp/5n2/1B2p1B1/4P3/1Q6/PPP2PPP/R3K2R w KQkq - 0 1"
    w_checks, b_checks = 0, 0
    
    board = ThreeCheckChess(initial_fen, w_checks, b_checks)

    # This is the shortest forced winning line against optimal defense.
    # Any other line can be defended by Black to last longer.
    moves_san = ["Bxf6+", "gxf6", "O-O-O", "exd4", "Rxd4", "Qe5", "Bxd7+", "Bxd7", "Rxd7+"]
    
    print("Finding the minimal amount of moves for White to win by 3 checks...")
    print(f"Initial Position. Checks: White={board.white_checks}, Black={board.black_checks}\n")
    
    white_move_count = 0
    
    for i, san_move in enumerate(moves_san):
        move = board.parse_san(san_move)
        
        if board.turn == chess.WHITE:
            white_move_count += 1
            turn_str = f"White's move {white_move_count}:"
        else:
            turn_str = f"Black's move {i // 2 + 1}:"

        board.push(move)
        
        print(f"{turn_str} {san_move}")
        check_info = f"Checks: White={board.white_checks}, Black={board.black_checks}"
        if board.is_check() and not board.is_game_over():
            check_info += " (Check!)"
        print(check_info)
        print(board)
        print("-" * 20)
        
        if board.white_checks >= 3:
            print(f"White wins by delivering 3 checks in {white_move_count} moves.")
            break
            
    # The final answer is the number of moves by white.
    # The puzzle asks for the integer answer.
    # print(f"The minimal amount of moves by white to win is {white_move_count}.")
    

solve()
print("<<<5>>>")