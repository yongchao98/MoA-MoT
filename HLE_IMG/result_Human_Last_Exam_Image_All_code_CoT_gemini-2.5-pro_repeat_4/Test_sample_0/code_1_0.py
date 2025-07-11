import sys

def solve_chess_puzzle():
    """
    This function prints the solution to the chess puzzle.
    The puzzle is a mate in 2 for Black, without moving the queens.

    The key move is Ng3+. This forces White into one of two replies,
    both of which lead to an immediate checkmate on Black's next move.

    Variation 1:
    1. ... Ng3+
    2. Kh1 ... (White's king moves to h1)
    3. ... Rg8# (Black's rook delivers mate)

    Variation 2:
    1. ... Ng3+
    2. hxg3 ... (White's pawn captures the knight)
    3. ... Rf1# (Black's rook delivers mate, which is possible because White's bishop on c1 is trapped and cannot capture the rook)

    The puzzle asks for the sequence of black's moves.
    """
    first_move = "Ng3+"
    second_move_variation_1 = "Rg8#"
    second_move_variation_2 = "Rf1#"

    # The problem asks for the sequence. We present the first move,
    # and then the two possible mating second moves.
    # The notation "move1 / move2" is used to show alternatives.
    print(f"1... {first_move}  2... {second_move_variation_1} / {second_move_variation_2}")

solve_chess_puzzle()