def solve_go_problem():
    """
    This function identifies the two optimal moves for Black and prints them.

    Analysis:
    1.  An endangered Black group is at (1,3) and (2,3). It needs to survive.
    2.  A White group of two stones at (3,3) and (4,3) is in atari (one liberty) at (4,2).
    3.  Move A: Black plays at (4,2). This immediately captures the two White stones.
        This helps the Black group survive by increasing its liberties from 2 to 3.
    4.  Move B: Black plays at (2,4). This connects the endangered Black group to a safe group,
        ensuring its survival. This also guarantees the capture of the White stones on the
        next turn, as they remain in atari and cannot be saved by White.
    5.  Both moves satisfy the problem's conditions. The answers are (2,4) and (4,2).
    6.  The answers must be sorted lexicographically.
    """
    # The two correct moves as (row, col) tuples
    answers = [(2, 4), (4, 2)]

    # Sort the answers lexicographically
    answers.sort()

    # Format the output string as requested
    # e.g., (r1, c1), (r2, c2)
    output_string = f"({answers[0][0]}, {answers[0][1]}), ({answers[1][0]}, {answers[1][1]})"
    
    print(output_string)

solve_go_problem()