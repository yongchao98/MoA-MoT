def solve_chess_puzzle():
    """
    Analyzes the chess position at move 130 from Carlsen-Nepomniachtchi, Game 6, 2021,
    and identifies the drawing move for Black.
    """

    explanation = """
    Step-by-step analysis to find the drawing move for Black:

    1.  **Position Analysis (before Black's 130th move):**
        White's pieces (King on h3, Rook on f5, Knight on h5) and the advanced pawn on e5 create a winning attack against Black's king on e8. Black's queen on a2 is too far away to defend conventionally. Black's only chance for survival is to create immediate counterplay with a perpetual check.

    2.  **The Losing Blunder:**
        In the actual game, Black played `130... Qe6`. This move fails because it allows White's king to improve its position with `131. Kh4`, after which White's attack is decisive, and Black cannot prevent the loss.

    3.  **The Drawing Combination:**
        The correct path to a draw is to initiate a perpetual check. Several queen checks can achieve this, but let's analyze option C, `130... Qg2`.

        This move forces the following sequence:

        -   Move 130 for Black: **Qg2+**
            This puts the white king in check.
        -   Move 131 for White: **Kh4**
            This is the only legal move for White's king.
        -   Move 131 for Black: **Qh2+**
            Black continues the series of checks.
        -   Move 132 for White: **Kg4**
            Again, this is the only legal move for White's king.
        -   Move 132 for Black: **Qg2+**
            Black forces the king back to h4, and the sequence repeats.

    4.  **The "Final Equation" of Moves:**
        The drawing sequence demonstrates a forced repetition of moves from which White cannot escape. The sequence is:

        130... Qg2+
        131. Kh4
        131... Qh2+
        132. Kg4
        132... Qg2+
        ... leading to a draw.

    5.  **Conclusion:**
        Move `130... Qg2` successfully forces a draw by perpetual check. Other checks like `Qh2+` or `Qf2+` would also lead to a similar drawing pattern.
    """
    print(explanation)
    # The final answer is the move that initiates this drawing sequence.
    final_answer = "C"
    print(f"<<<{final_answer}>>>")

solve_chess_puzzle()