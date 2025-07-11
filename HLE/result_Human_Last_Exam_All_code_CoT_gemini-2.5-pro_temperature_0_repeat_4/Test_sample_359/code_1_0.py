import sys

def solve_koth_puzzle():
    """
    Analyzes the King of the Hill chess position and determines the minimum number of moves for White to win.
    
    FEN: 8/2k5/5pn1/1Pp1pNpp/3PP3/4K1B1/8/8 w - - 0 43
    Winning condition: Move king to d4, e4, d5, or e5.
    """

    # Step 1: Analyze the initial position.
    # White King: e3
    # White Pawns: d4, e4
    # Black King: c7
    # Black Pawn: c5
    # White Bishop: g3
    # The central squares d4 and e4 are blocked by White's own pawns.
    # White needs to clear a path for the king to the center.

    # Step 2: Identify White's optimal move sequence.
    # The most powerful move for White is dxc5.
    
    print("The winning plan for White involves a 2-move sequence.")
    print("-" * 30)

    # Move 1 for White
    move_number_1 = 1
    white_move_1 = "dxc5+"
    print(f"Move {move_number_1} (White): {white_move_1}")
    print("The White pawn on d4 captures the Black pawn on c5.")
    print("This action has two critical effects:")
    print("1. It vacates the central square d4.")
    print("2. It is a discovered check (+) from the White Bishop on g3 to the Black King on c7.")
    
    # Move 1 for Black
    black_move_1 = "...Kb8"
    print(f"\nMove {move_number_1} (Black): {black_move_1} (or any other legal king move)")
    print("Black is in check and is forced to move the king. Black cannot block the check or capture the bishop.")
    print("The specific square the king moves to does not prevent White's next move.")

    # Move 2 for White
    move_number_2 = 2
    white_move_2 = "Kd4"
    print(f"\nMove {move_number_2} (White): {white_move_2}")
    print("The d4 square is now empty and is not attacked by any Black piece.")
    print("White moves the King from e3 to d4.")
    print("By moving the King to a central square, White wins the game.")
    print("-" * 30)

    # Final Conclusion
    total_moves_for_win = 2
    print(f"\nConclusion: White can force a win in {total_moves_for_win} moves.")

solve_koth_puzzle()