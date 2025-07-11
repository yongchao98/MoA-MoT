import collections

def solve_chess_puzzle():
    """
    This script determines the minimum piece configuration for a specific
    Diagonal Corridor Mate, then prints the reasoning and the final answer.
    """
    # Step 1: Define the problem parameters and piece values
    print("Goal: Find the minimal piece configuration for a Diagonal Corridor Mate.")
    print("Kings: White on a1, Black on h8. White to deliver mate.")
    print("Constraints:")
    print("1. Use the minimum number of additional pieces.")
    print("2. If tied, use the minimum total value of the additional pieces.")
    print("-" * 30)

    piece_values = {'Pawn': 1, 'Knight': 3, 'Bishop': 3, 'Rook': 5, 'Queen': 9}

    # Step 2: Analyze the checkmate requirements
    print("Analysis of the Mate:")
    print("The Black King on h8 has three escape squares: g7, g8, and h7.")
    print("For a checkmate, these squares must be blocked or attacked.")
    print("The check must come from a piece on the a1-h8 diagonal.")
    print("-" * 30)

    # Step 3: Determine the most efficient piece configuration
    print("Finding the Minimal Solution:")
    print("We seek a 3-piece solution (1 White attacker, 2 Black blockers), as this is the theoretical minimum.")
    print("The two Black pieces will be Pawns to minimize value, blocking h7 and g8.")
    print("The single White piece must check h8 AND control the last escape square, g7.")
    print("\nEvaluating potential White pieces:")

    options = []
    pawn_value = piece_values['Pawn']

    # Option A: Using a White Bishop
    bishop_value = piece_values['Bishop']
    option_A_value = bishop_value + pawn_value + pawn_value
    options.append({
        'name': 'Bishop Solution',
        'white_pieces': ['White Bishop'],
        'black_pieces': ['Black Pawn', 'Black Pawn'],
        'total_value': option_A_value
    })
    print(f"- Using a White Bishop: It can attack h8 and g7 from the same diagonal.")
    print(f"  Equation for this solution's value: {bishop_value} (Bishop) + {pawn_value} (Pawn) + {pawn_value} (Pawn) = {option_A_value}")

    # Option B: Using a White Queen
    queen_value = piece_values['Queen']
    option_B_value = queen_value + pawn_value + pawn_value
    options.append({
        'name': 'Queen Solution',
        'white_pieces': ['White Queen'],
        'black_pieces': ['Black Pawn', 'Black Pawn'],
        'total_value': option_B_value
    })
    print(f"- Using a White Queen: It can also attack both h8 and g7.")
    print(f"  Equation for this solution's value: {queen_value} (Queen) + {pawn_value} (Pawn) + {pawn_value} (Pawn) = {option_B_value}")
    print("-" * 30)

    # Step 4: Select the best option based on the criteria
    print("Selecting the Optimal Solution:")
    best_solution = min(options, key=lambda x: x['total_value'])
    print(f"Both solutions use 3 pieces. The tie-breaker is minimum piece value.")
    print(f"The '{best_solution['name']}' is optimal with a total value of {best_solution['total_value']}.")
    print("-" * 30)

    # Step 5: Format and print the final answer string
    print("Final Answer Construction:")
    white_pieces_list = sorted(best_solution['white_pieces'])
    black_pieces_list = sorted(best_solution['black_pieces'])
    
    # Joining the lists with a comma
    final_answer_list = white_pieces_list + black_pieces_list
    final_answer_str = ", ".join(final_answer_list)
    
    print("The required list of additional pieces is:")
    print(final_answer_str)

solve_chess_puzzle()