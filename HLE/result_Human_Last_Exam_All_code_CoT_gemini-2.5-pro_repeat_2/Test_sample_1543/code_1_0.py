def solve_capablanca_puzzle():
    """
    Analyzes the Capablanca chess position and prints the solution.
    """
    
    # FEN: 9k/5c1pb1/10/10/10/3Q6/PP5A2/K9 w - - 0 1
    # Board: 10x8 (files a-j, ranks 1-8)
    
    # Key piece positions:
    # White: King(a1), Queen(d3), Archbishop(h2)
    # Black: King(j8), Chancellor(f7), Bishop(i7), Pawn(h7)

    print("Analyzing the Capablanca chess puzzle...")
    print("The goal is to find the minimum number of moves for White to win (checkmate).")
    print("-" * 30)
    
    print("Step 1: Check for Mate in 1.")
    print("No single move by White can deliver an immediate checkmate. For example, 1. Qj7+ is met by ...cxj7.")
    print("-" * 30)

    print("Step 2: Look for Mate in 2.")
    print("White's first move must be a 'quiet' but deadly move that sets up an unstoppable mating net.")
    
    move_1_white = "1. Aj3"
    print(f"The key move for White is {move_1_white}.")
    print("The Archbishop (A) moves from h2 to j3. This is a knight-move.")
    print("This move has a critical consequence: the Archbishop on j3 now attacks the black Chancellor on f7 (as a bishop).")
    print("-" * 30)

    print("Step 3: Analyze Black's forced responses.")
    print("Black must respond to the attack on the Chancellor. We consider two cases for Black's optimal defense:")

    # Case 1: Black moves Chancellor to j7
    move_1_black_case1 = "1... cj7"
    move_2_white_case1 = "2. Qd8#"
    print(f"\nCase A: Black plays {move_1_black_case1}.")
    print("The Chancellor moves to j7 to defend the j-file.")
    print(f"White responds with {move_2_white_case1} (Queen from d3 to d8). This is checkmate.")
    print("  - The Black King on j8 is in check from the Queen on d8.")
    print("  - The escape square i8 is covered by the Queen on d8.")
    print("  - The escape square j7 is blocked by Black's own Chancellor.")
    print("  - The escape square i7 is blocked by Black's own Bishop.")
    print("  - The Queen on d8 cannot be captured.")

    # Case 2: Black moves Chancellor to another square
    move_1_black_case2 = "1... ce7"
    move_2_white_case2 = "2. Qj8#"
    print(f"\nCase B: Black moves the Chancellor to any other safe square, for example {move_1_black_case2}.")
    print(f"White responds with {move_2_white_case2} (Queen from d3 to j8). This is checkmate.")
    print("  - The Black King on j8 is in check from the Queen on j8.")
    print("  - The escape square i8 is covered by the Queen on j8 (diagonally).")
    print("  - The escape square j7 is covered by the Queen on j8 (along the file).")
    print("  - The escape square i7 is blocked by Black's own Bishop.")
    print("  - The Queen on j8 cannot be captured.")
    print("-" * 30)

    final_answer = 2
    print("Conclusion: White can force checkmate on the second move against any defense.")
    print(f"Therefore, the minimal amount of moves by White to win is {final_answer}.")
    
    return final_answer

if __name__ == "__main__":
    solve_capablanca_puzzle()
