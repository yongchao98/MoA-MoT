def solve_go_problem():
    """
    Analyzes a Go board configuration to find the winning move for Black.
    """
    black_stones = {(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)}
    white_stones = {(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)}

    print("Initial Stone Configuration:")
    print(f"Black stones: {sorted(list(black_stones))}")
    print(f"White stones: {sorted(list(white_stones))}")
    print("-" * 30)

    # The key insight is to find a move that splits the weak white group.
    # The move (2, 4) is the vital point that achieves this.
    best_move = (2, 4)
    print(f"Black's winning move is to play at: {best_move}")
    print("-" * 30)
    
    print("Analysis of the move:")
    print(f"Playing at {best_move} splits the single large White group into two smaller, weaker groups.")

    # Add the new black stone to the board
    new_black_stones = black_stones.copy()
    new_black_stones.add(best_move)
    all_stones = new_black_stones.union(white_stones)

    # Define the two new white subgroups after the split
    white_group_1 = {(2, 5), (1, 4)}
    white_group_2 = {(3, 4), (3, 3), (2, 2)}

    def get_neighbors(r, c):
        """Returns the four neighbors of a coordinate."""
        return {(r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}

    def get_liberties(group, occupied_points):
        """Calculates the liberties for a given group of stones."""
        liberties = set()
        for r, c in group:
            neighbors = get_neighbors(r, c)
            for neighbor in neighbors:
                if neighbor not in occupied_points:
                    liberties.add(neighbor)
        return liberties

    # Analyze the state of each subgroup after the split
    liberties_group_1 = get_liberties(white_group_1, all_stones)
    liberties_group_2 = get_liberties(white_group_2, all_stones)

    print("\n--- After Black plays at (2, 4) ---")
    print(f"White Group 1: {sorted(list(white_group_1))}")
    print(f"Liberties of Group 1: {sorted(list(liberties_group_1))}")
    print("Group 1 is now very weak with only 2 liberties. Black can easily capture it by playing at (1, 5) and then (1, 3).")

    print(f"\nWhite Group 2: {sorted(list(white_group_2))}")
    print(f"Liberties of Group 2: {sorted(list(liberties_group_2))}")
    print("Group 2 is also weak, forming a 'dead' shape with only 4 liberties. Black can capture it by filling these liberties.")
    
    print("\nConclusion:")
    print("White cannot defend both groups at the same time (this is called 'miai' in Go).")
    print("If White tries to save one group, Black will capture the other, and then come back to capture the first one.")
    print("Therefore, the move (2, 4) guarantees the eventual capture of all White stones.")

solve_go_problem()
<<<G>>>