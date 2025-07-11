import sys

def solve_chess_puzzle():
    """
    This function prints the step-by-step move sequence for the forced checkmate.
    """
    # The forced mating sequence involves a series of checks leading to a smothered mate.
    # The sequence is 1... Nf2+ 2. Kg1 Nh3+ 3. Kh1 Qg1+ 4. Rxg1 Nf2#

    # Define the moves in the sequence
    moves = [
        ("1...", "Nf2+"),
        ("2.", "Kg1", "Nh3+"),
        ("3.", "Kh1", "Qg1+"),
        ("4.", "Rxg1", "Nf2#")
    ]
    
    # Print the moves step by step
    print("The best sequence of moves is:")
    
    final_equation_parts = []
    
    # 1. Black's first move
    move_number_black = moves[0][0]
    black_move = moves[0][1]
    line = f"{move_number_black} {black_move}"
    final_equation_parts.append(line)
    
    # 2. White's reply and Black's second move
    move_number_white = moves[1][0]
    white_move = moves[1][1]
    black_move = moves[1][2]
    line = f"{move_number_white} {white_move} {black_move}"
    final_equation_parts.append(line)
    
    # 3. White's reply and Black's third move
    move_number_white = moves[2][0]
    white_move = moves[2][1]
    black_move = moves[2][2]
    line = f"{move_number_white} {white_move} {black_move}"
    final_equation_parts.append(line)
    
    # 4. White's reply and Black's final mating move
    move_number_white = moves[3][0]
    white_move = moves[3][1]
    black_move = moves[3][2]
    line = f"{move_number_white} {white_move} {black_move}"
    final_equation_parts.append(line)

    final_equation = " ".join(final_equation_parts)
    print("\nThe final equation representing the sequence is:")
    print(final_equation)


solve_chess_puzzle()
# The final response should have the answer tag at the end.
# I will check if my answer matches the options provided.
# A. 1... Qxe1+, 2. Rxe1 Nf2+, 3. Kg1 Nh3# (Incorrect line)
# B. 1... Ng3+ 2. hxg3 Qxg3 3. Qe1 Qh4# (Not forced)
# C. 1... Qg1+ 2. Rxg1 Nf2# (Shorter, but the chosen line is also a valid pattern)
# D. 1...Qg1+ 2.Kxg1 Ne2+ 3.Kf1 N4g3+ 4.hxg3# (Illegal move Kxg1)
# E. 1... Nf2+ 2.Kg1 Nh3+ 3.Kh1 Qg1+ 4.Rxg1 Nf2# (This is the sequence I've described)
# My solution corresponds to option E.
sys.stdout.write("\n<<<E>>>")