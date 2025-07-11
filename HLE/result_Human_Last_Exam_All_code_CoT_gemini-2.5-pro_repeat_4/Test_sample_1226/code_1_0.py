def solve_chess_puzzle():
    """
    Calculates the maximum material points for a mate in not less than 6 moves.
    """
    piece_values = {
        'Queen': 9,
        'Rook': 5,
        'Bishop': 3,
        'Knight': 3,
        'Pawn': 1
    }

    # 1. Material for the core mating mechanism.
    # A Queen at a3 and King at c3 force the bK to shuttle between a1 and b1.
    shuttle_queen_count = 1
    shuttle_queen_val = shuttle_queen_count * piece_values['Queen']

    # A mating Queen starts at h8 to deliver mate at h2 on move 6.
    mating_queen_count = 1
    mating_queen_val = mating_queen_count * piece_values['Queen']

    # 5 Queens block the mating path on h3, h4, h5, h6, h7.
    # White must spend 5 moves to move them.
    num_blocking_queens = 5
    blocking_queens_val = num_blocking_queens * piece_values['Queen']

    # 2. Calculate the number of available squares for "filler" pieces.
    total_squares = 64
    
    # Squares used by the solution mechanism (initial positions, paths, destinations).
    # Initial: bK@a1, wK@c3, wQ@a3, wQ@h8, wQ@h7, wQ@h6, wQ@h5, wQ@h4, wQ@h3 (9 squares)
    # Paths/Dests: bK path to b1, wQ mate on h2, plus 5 safe squares for blockers (7 squares)
    solution_squares_count = 16

    # Squares from which a Queen could attack the bK shuttle squares (a1 or b1),
    # creating a shorter mate or an illegal position.
    # Files 'a'/'b'(16) + Rank '1'(6) + Diagonals from a1/b1 (11) = 33 squares.
    attacking_squares_count = 33
    
    # The sets of solution squares and attacking squares overlap on 4 squares (a1, b1, a3, c3).
    overlap_count = 4
    
    unavailable_squares = solution_squares_count + attacking_squares_count - overlap_count
    available_filler_squares = total_squares - unavailable_squares
    
    # 3. Calculate points from filler Queens on available squares.
    num_filler_queens = available_filler_squares
    filler_queens_val = num_filler_queens * piece_values['Queen']

    # 4. Sum all the points.
    total_value = shuttle_queen_val + mating_queen_val + blocking_queens_val + filler_queens_val

    # Print the breakdown of the calculation.
    print("The maximum point value is derived from the following setup:")
    print(f"- Shuttle Queen: {shuttle_queen_val} points")
    print(f"- Mating Queen: {mating_queen_val} points")
    print(f"- Blocking Queens: {blocking_queens_val} points")
    print(f"- Filler Queens: {filler_queens_val} points ({num_filler_queens} Queens on safe squares)")
    print("\nFinal Point Calculation:")
    print(f"{shuttle_queen_val} + {mating_queen_val} + {blocking_queens_val} + {filler_queens_val} = {total_value}")

solve_chess_puzzle()