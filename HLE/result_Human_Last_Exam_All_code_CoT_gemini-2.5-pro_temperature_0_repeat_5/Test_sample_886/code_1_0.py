def solve_chess_puzzle():
    """
    Calculates the number of distinct checkmate positions achievable by a single
    super-piece against a lone king without assistance.
    """
    print("Step 1: Analyze the conditions for an unassisted checkmate.")
    print("An unassisted checkmate by a single super-piece against a lone king is only possible if the enemy king is in a corner.")
    print("A chessboard has 4 corners.")

    print("\nStep 2: Identify the unique mating patterns based on the super-pieces' moves.")
    print("We analyze the unique squares from which a super-piece can deliver a corner checkmate.")
    
    # Mating patterns for a king in a corner (e.g., h8):
    # 1. A piece with Queen moves (Queen, Amazon, Q+K) can mate from 2 squares (e.g., f7, g6).
    queen_like_mates = 2
    print(f"Pattern A (Queen-like): A piece with Queen moves can mate from {queen_like_mates} distinct squares.")
    
    # 2. A piece with Bishop+Knight moves (Archbishop, Amazon) can mate from 1 additional square (e.g., f6).
    archbishop_like_mates = 1
    print(f"Pattern B (Archbishop-like): A piece with Bishop+Knight moves adds {archbishop_like_mates} new mating square.")

    # 3. A piece with Pawn+Rook/King/Queen moves can mate from 1 square (e.g., g7), but this is conditional.
    pawn_like_mates = 1
    print(f"Pattern C (Pawn-synergy): A piece with Pawn moves adds {pawn_like_mates} new mating square, but only if its pawn-attack can target the corner.")

    print("\nStep 3: Calculate mating positions for a White super-piece vs. a Black King.")
    print("A White piece's pawn-move can only attack a Black King on the 8th rank (corners a8, h8).")
    
    corners_pawn_possible = 2
    mates_per_pawn_corner = queen_like_mates + archbishop_like_mates + pawn_like_mates
    print(f"For the {corners_pawn_possible} corners where a pawn-like mate is possible, there are {queen_like_mates} + {archbishop_like_mates} + {pawn_like_mates} = {mates_per_pawn_corner} mating squares each.")

    print("A White piece's pawn-move cannot attack a Black King on the 1st rank (corners a1, h1).")
    corners_pawn_impossible = 2
    mates_per_no_pawn_corner = queen_like_mates + archbishop_like_mates
    print(f"For the {corners_pawn_impossible} corners where a pawn-like mate is not possible, there are {queen_like_mates} + {archbishop_like_mates} = {mates_per_no_pawn_corner} mating squares each.")

    white_mates_black_total = (corners_pawn_possible * mates_per_pawn_corner) + (corners_pawn_impossible * mates_per_no_pawn_corner)
    print(f"Total positions for White mating Black: ({corners_pawn_possible} * {mates_per_pawn_corner}) + ({corners_pawn_impossible} * {mates_per_no_pawn_corner}) = {white_mates_black_total}")

    print("\nStep 4: Calculate mating positions for a Black super-piece vs. a White King.")
    print("By symmetry, the number of positions is the same. A Black piece can perform a pawn-like mate on the 1st rank.")
    black_mates_white_total = white_mates_black_total
    print(f"Total positions for Black mating White: {black_mates_white_total}")

    print("\nStep 5: Calculate the total number of distinct checkmate positions.")
    total_positions = white_mates_black_total + black_mates_white_total
    print("The total is the sum of positions for both scenarios (White mates Black, and Black mates White).")
    
    print("\nFinal Calculation:")
    print(f"{white_mates_black_total} + {black_mates_white_total} = {total_positions}")

solve_chess_puzzle()
<<<28>>>