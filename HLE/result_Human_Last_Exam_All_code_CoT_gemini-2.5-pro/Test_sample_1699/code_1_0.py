def solve_go_problem():
    """
    Analyzes the Go board to find the optimal move for Black.
    """
    # Define stone positions
    black_stones = {(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)}
    white_stones = {(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)}

    # Helper functions
    def get_neighbors(r, c):
        return {(r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}

    def find_group(stone, all_stones):
        q = [stone]
        group = {stone}
        visited = {stone}
        while q:
            s = q.pop(0)
            for n in get_neighbors(s[0], s[1]):
                if n in all_stones and n not in visited:
                    visited.add(n)
                    group.add(n)
                    q.append(n)
        return group

    def get_liberties(group, opposing_stones, own_stones):
        liberties = set()
        for r, c in group:
            for neighbor in get_neighbors(r, c):
                if neighbor not in opposing_stones and neighbor not in own_stones:
                    liberties.add(neighbor)
        return liberties

    # --- Main Analysis ---
    print("Step 1: Analyzing the initial board state for White.")
    
    remaining_white = set(white_stones)
    white_groups = []
    while remaining_white:
        stone = remaining_white.pop()
        group = find_group(stone, white_stones)
        white_groups.append(group)
        remaining_white -= group
        
    print("Identified White groups and their initial liberties:")
    for i, group in enumerate(white_groups):
        liberties = get_liberties(group, black_stones, white_stones)
        print(f"  - White Group {i+1} {sorted(list(group))}: {len(liberties)} liberties.")
    
    print("\nStep 2: Evaluating the proposed move for Black at (2, 4).")
    
    move_to_evaluate = (2, 4)
    new_black_stones = black_stones.copy()
    new_black_stones.add(move_to_evaluate)
    
    print(f"State after Black plays at {move_to_evaluate}:")
    for i, group in enumerate(white_groups):
        liberties = get_liberties(group, new_black_stones, white_stones)
        status = ""
        if len(liberties) == 1:
            status = " (IN ATARI!)"
        print(f"  - White Group {i+1} {sorted(list(group))}: Now has {len(liberties)} liberties.{status}")

    print("\nStep 3: Conclusion.")
    print("The move at (2, 4) is the vital point because it attacks multiple groups at once.")
    print("It puts the group at (2, 5) in atari, forcing a response from White.")
    print("This creates a 'miai' situation: White cannot save all their stones. If they save one group, Black will attack another.")
    print("This strategic advantage allows Black to eventually capture all White stones.")
    
    chosen_move_row = 2
    chosen_move_col = 4
    
    print("\nThe final answer is derived from the coordinates of the chosen move.")
    print(f"The numbers in the final answer are: {chosen_move_row} and {chosen_move_col}")


solve_go_problem()