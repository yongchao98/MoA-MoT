import collections

def solve_chess_puzzle():
    """
    This function solves the Diagonal Corridor Mate puzzle based on the given constraints.
    It determines the minimum number of pieces required and then considers the minimum piece value.
    """

    # Piece values for sorting and evaluation
    piece_values = {
        "Pawn": 1,
        "Knight": 3,
        "Bishop": 3,
        "Rook": 5,
        "Queen": 9,
        "King": 1000  # Not used in calculation for additional pieces
    }

    # Step 1: Define the core of the checkmate.
    # The Black King is on h8. White mates along the a1-h8 diagonal.
    # The most efficient mating piece for White is a Bishop.
    white_pieces = ["Bishop"]
    white_value = piece_values["Bishop"]

    # Step 2: Block the Black King's escape squares (g7, h7, g8).
    # We need to find the minimal set of Black pieces to do this.
    # A single Black Queen on g7 can block all three escape squares.
    # This leads to a solution with only 2 total additional pieces.
    
    # Candidate Solution 1 (2 pieces):
    solution_1_black_pieces = ["Queen"]
    solution_1_black_value = piece_values["Queen"]
    solution_1_total_pieces = len(white_pieces) + len(solution_1_black_pieces)
    solution_1_is_valid = solution_1_black_value > white_value
    
    # Let's consider an alternative with more pieces to prove the 2-piece solution is minimal.
    # For example, a Black Rook on g7 (blocks g7, g8) and a Black Pawn on h7 (blocks h7).
    
    # Candidate Solution 2 (3 pieces):
    solution_2_black_pieces = ["Rook", "Pawn"]
    solution_2_black_value = piece_values["Rook"] + piece_values["Pawn"]
    solution_2_total_pieces = len(white_pieces) + len(solution_2_black_pieces)
    solution_2_is_valid = solution_2_black_value > white_value

    # Step 3: Compare solutions.
    # The primary condition is the minimum number of pieces.
    # Solution 1 has 2 pieces. Solution 2 has 3 pieces.
    # Therefore, Solution 1 is the correct one.
    
    final_white_pieces = white_pieces
    final_black_pieces = solution_1_black_pieces
    
    # Step 4: Format the output as requested.
    # Sort each list by piece value, then by name alphabetically as a tie-breaker.
    final_white_pieces.sort(key=lambda p: (piece_values[p], p))
    final_black_pieces.sort(key=lambda p: (piece_values[p], p))

    # Combine into a single comma-separated string.
    white_pieces_str = [f"White {p}" for p in final_white_pieces]
    black_pieces_str = [f"Black {p}" for p in final_black_pieces]
    
    final_answer = ", ".join(white_pieces_str + black_pieces_str)

    print("To construct the Diagonal Corridor Mate with the minimum number of pieces:")
    print("1. Mating Piece (White): A Bishop is needed to attack the Black King on h8 along the long diagonal. This is the lowest value piece that can perform this role.")
    print("2. Blocking Pieces (Black): The King's escape squares (g7, h7, g8) must be blocked. A single Black Queen, placed on g7, can block all three squares simultaneously.")
    print("3. Evaluation: This setup uses 2 additional pieces (1 White, 1 Black). Black's piece value (Queen=9) is greater than White's (Bishop=3), satisfying all conditions.")
    print("\nFinal list of additional pieces:")
    print(final_answer)

solve_chess_puzzle()
<<<White Bishop, Black Queen>>>