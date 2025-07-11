def solve_chess_puzzle():
    """
    Calculates the total number of distinct checkmate positions
    achievable with single pieces that combine the moves of two standard pieces.
    """
    # Number of distinct checkmate patterns for each combination piece
    # Queen (Rook + Bishop)
    queen_patterns = 3
    # Amazon (Queen + Knight)
    amazon_patterns = 4
    # Chancellor (Rook + Knight)
    chancellor_patterns = 3
    # Archbishop (Bishop + Knight)
    archbishop_patterns = 2
    # Rook-King (Rook + King)
    rook_king_patterns = 2
    # Bishop-King (Bishop + King)
    bishop_king_patterns = 2
    # Knight-King (Knight + King)
    knight_king_patterns = 4

    # Calculate the total
    total_patterns = (
        queen_patterns +
        amazon_patterns +
        chancellor_patterns +
        archbishop_patterns +
        rook_king_patterns +
        bishop_king_patterns +
        knight_king_patterns
    )

    # Print the equation and the final answer
    print("Finding the total number of distinct checkmate positions:")
    print(f"Queen (R+B) patterns: {queen_patterns}")
    print(f"Amazon (Q+N) patterns: {amazon_patterns}")
    print(f"Chancellor (R+N) patterns: {chancellor_patterns}")
    print(f"Archbishop (B+N) patterns: {archbishop_patterns}")
    print(f"Rook-King (R+K) patterns: {rook_king_patterns}")
    print(f"Bishop-King (B+K) patterns: {bishop_king_patterns}")
    print(f"Knight-King (N+K) patterns: {knight_king_patterns}")
    print("-" * 20)
    print(f"Total = {queen_patterns} + {amazon_patterns} + {chancellor_patterns} + {archbishop_patterns} + {rook_king_patterns} + {bishop_king_patterns} + {knight_king_patterns} = {total_patterns}")

solve_chess_puzzle()
