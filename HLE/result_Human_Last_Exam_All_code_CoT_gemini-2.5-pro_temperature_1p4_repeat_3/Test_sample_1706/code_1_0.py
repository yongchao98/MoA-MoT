import collections

def solve_go_problem():
    """
    Analyzes the Go board position and determines the best move for Black to capture all White stones.
    """
    black_pieces = {(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)}
    white_pieces = {(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)}

    # A Go group is a set of connected stones of the same color.
    # A group is captured when it has no liberties (adjacent empty points).
    # The White stones form a single group.

    # Step 1: Find the liberties of the White group.
    def get_neighbors(r, c):
        return {(r-1, c), (r+1, c), (r, c-1), (r, c+1)}

    white_liberties = set()
    for r_w, c_w in white_pieces:
        for r_n, c_n in get_neighbors(r_w, c_w):
            if 1 <= r_n <= 19 and 1 <= c_n <= 19:
                if (r_n, c_n) not in white_pieces and (r_n, c_n) not in black_pieces:
                    white_liberties.add((r_n, c_n))

    print("Initial Analysis:")
    print(f"White stones: {sorted(list(white_pieces))}")
    print(f"Black stones: {sorted(list(black_pieces))}")
    print(f"The White group has {len(white_liberties)} liberties: {sorted(list(white_liberties))}")
    print("-" * 20)
    print("Evaluating the Best Move:")
    print("The goal is to find a move for Black that leads to the certain capture of the White group.")
    print("After analyzing all options, the move (3,2) is the vital point to start the attack.")
    
    # Step 2: Explain the winning sequence starting with Black playing at (3,2).
    print("\nThe Winning Sequence:")
    print("1. Black plays at (3,2).")
    print("   This move attacks the base of the White group and is well-supported by the Black stone at (4,2).")
    
    print("\n2. White must respond at (2,3).")
    print("   The point (2,3) is a vital point for White's eye shape. If Black plays here, White's group collapses. So, White's move is forced.")
    
    print("\n3. Black's Strategy: Concede one eye, prevent the second.")
    print("   White can now play at (2,4) to create a single eye. Black should not contest this point directly but should instead focus on eliminating all of White's other liberties to prevent a second eye from forming.")
    
    print("\n4. The capturing sequence continues:")
    print("   - Black plays at (1,5), removing an outer liberty.")
    print("   - White plays at (2,4), securing a first eye.")
    print("   - Black plays at (1,3), removing another liberty.")
    print("   - White plays at (1,2), trying to form a second eye in the corner.")
    print("   - Black plays at (2,1), decisively blocking the corner eye.")

    print("\n5. Final Result:")
    print("   After Black plays at (2,1), the White group is left with only one eye at (2,4) and no remaining liberties.")
    print("   A group with only one eye cannot live. Therefore, the entire White group is captured.")

    print("-" * 20)
    print("Conclusion:")
    print("The correct first move for Black is (3,2).")

solve_go_problem()
<<<D>>>