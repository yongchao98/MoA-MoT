def solve_chess_material_puzzle():
    """
    Calculates the maximum score of white material for a mate in exactly 6 moves.
    This function codifies the logic developed above.
    """
    # Piece values
    QUEEN = 9
    ROOK = 5
    BISHOP = 3
    PAWN = 1

    # 1. Score from the core trapping pieces
    # Black King at a1, White King at c2.
    # Trap pieces: Pawns at a2, b2, and a Bishop at b3.
    pawn_trapper1_score = PAWN
    pawn_trapper2_score = PAWN
    bishop_trapper_score = BISHOP
    trapping_score = pawn_trapper1_score + pawn_trapper2_score + bishop_trapper_score
    print("Score of trapping pieces (2 Pawns, 1 Bishop): {}".format(trapping_score))

    # 2. Score from the core mating apparatus
    # Mating Rook at h8, blocked by 5 Queens on b8, c8, d8, e8, f8.
    mating_rook_score = ROOK
    blocking_queens_score = QUEEN * 5
    mating_apparatus_score = mating_rook_score + blocking_queens_score
    print("Score of mating apparatus (1 Rook, 5 Queens): {}".format(mating_apparatus_score))
    
    core_score = trapping_score + mating_apparatus_score
    print("Total score of the 9 core White pieces: {}".format(core_score))
    
    # 3. Score from filling the rest of the board
    total_squares = 64
    # WK, BK, P, P, B, R, Q, Q, Q, Q, Q = 12 pieces
    core_pieces_count = 12
    empty_squares = total_squares - core_pieces_count
    
    # On squares forbidden for Queen/Bishop (c4, d5, e6, f7, g8), place Rooks.
    diag_forbidden_count = 5
    filler_rooks_score = diag_forbidden_count * ROOK
    print("\nScore of filler Rooks on diagonal: {}".format(filler_rooks_score))

    # On squares forbidden for Queen/Rook (b1, c1, d1, e1, f1, g1, h1), place Bishops.
    rank_forbidden_count = 7
    filler_bishops_score = rank_forbidden_count * BISHOP
    print("Score of filler Bishops on rank 1: {}".format(filler_bishops_score))

    # Fill all remaining squares with Queens.
    queen_filler_count = empty_squares - diag_forbidden_count - rank_forbidden_count
    filler_queens_score = queen_filler_count * QUEEN
    print("Number of filler Queens: {}".format(queen_filler_count))
    print("Score of filler Queens: {}".format(filler_queens_score))
    
    filler_score = filler_rooks_score + filler_bishops_score + filler_queens_score
    
    # 4. Final total score
    total_score = core_score + filler_score
    
    print("\nFinal Calculation:")
    print("Core piece score: {}".format(core_score))
    print("Filler piece score: {}".format(filler_score))
    print("The final equation is: {} (trapping) + {} (mating) + {} (filler rooks) + {} (filler bishops) + {} (filler queens)".format(
        trapping_score, mating_apparatus_score, filler_rooks_score, filler_bishops_score, filler_queens_score))
    print("\nTotal maximum material points:")
    print(total_score)

    return total_score

solve_chess_material_puzzle()
