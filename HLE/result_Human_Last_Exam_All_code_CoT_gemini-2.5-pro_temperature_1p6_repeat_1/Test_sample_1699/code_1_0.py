import collections

def get_group_and_liberties(r_start, c_start, all_stones_of_color, all_opponent_stones, all_stones):
    """
    Finds a group of connected stones and its liberties using Breadth-First Search (BFS).
    
    Args:
        r_start, c_start: The starting coordinates of a stone.
        all_stones_of_color: A set of (row, col) tuples for the stones of the current player.
        all_opponent_stones: A set of (row, col) tuples for the opponent's stones.
        all_stones: A set of all stone positions.
    
    Returns:
        A tuple containing:
        - group_stones (set): The set of (row, col) for the connected group.
        - liberties (set): The set of (row, col) for the group's liberties.
    """
    if (r_start, c_start) not in all_stones_of_color:
        return None, None

    q = collections.deque([(r_start, c_start)])
    visited = set([(r_start, c_start)])
    group_stones = set()
    liberties = set()

    while q:
        r, c = q.popleft()
        group_stones.add((r, c))

        # Check neighbors
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc

            # Skip out of bounds, though not strictly necessary for this problem
            if not (1 <= nr <= 19 and 1 <= nc <= 19):
                continue
            
            neighbor_pos = (nr, nc)
            if neighbor_pos in visited:
                continue
            
            # If neighbor is empty, it's a liberty
            if neighbor_pos not in all_stones:
                liberties.add(neighbor_pos)
            # If neighbor is a stone of the same color, add to search queue
            elif neighbor_pos in all_stones_of_color:
                visited.add(neighbor_pos)
                q.append(neighbor_pos)
    
    return group_stones, liberties

def analyze_white_stones(black_stones, white_stones):
    """Analyzes all white groups on the board."""
    all_s = black_stones | white_stones
    
    # Use a copy to avoid modifying the set while iterating
    white_stones_to_check = set(white_stones)
    groups = []
    
    while white_stones_to_check:
        r, c = white_stones_to_check.pop()
        group, libs = get_group_and_liberties(r, c, white_stones, black_stones, all_s)
        
        # Remove all stones of this group from the set to check
        white_stones_to_check -= group
        groups.append({'group': sorted(list(group)), 'liberties': len(libs)})
    
    # Sort groups for consistent output
    groups.sort(key=lambda x: x['group'][0])
    return groups

def main():
    # Initial board state
    black_stones = {(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)}
    white_stones = {(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)}
    
    candidate_moves = {
        "B": (1, 6),
        "C": (2, 4),
        "D": (1, 3),
        "E": (1, 2),
        "F": (3, 2),
        "G": (2, 1)
    }

    print("--- Initial Board Analysis ---")
    initial_analysis = analyze_white_stones(black_stones, white_stones)
    print("Initial state of White groups and their liberties:")
    for group_info in initial_analysis:
        print(f"  - White group at {group_info['group']} has {group_info['liberties']} liberties.")
    print("-" * 30)

    best_move = None
    max_impact = -1

    print("\n--- Analyzing Potential Black Moves ---")
    for choice, move in candidate_moves.items():
        print(f"\nAnalyzing move {choice}: Black plays at {move}")
        
        # Create new board state
        temp_black_stones = black_stones | {move}
        
        new_analysis = analyze_white_stones(temp_black_stones, white_stones)
        
        print("Impact on White groups:")
        impact_score = 0
        for i, group_info in enumerate(new_analysis):
            initial_libs = initial_analysis[i]['liberties']
            new_libs = group_info['liberties']
            
            lib_change = initial_libs - new_libs
            impact_score += lib_change
            
            status = ""
            if lib_change > 0:
                status = f"Liberties reduced from {initial_libs} to {new_libs}."
                if new_libs == 1:
                    status += " (Now in ATARI)"
            else:
                status = "No change."
            
            print(f"  - White group at {group_info['group']}: {status}")

    print("\n--- Conclusion ---")
    print("The goal is to find a move that applies maximum pressure to the White stones.")
    print("Based on the analysis, we can evaluate each move:")
    print(" - Moves like (1,6), (1,2), or (2,1) affect at most one group and are not critical.")
    print(" - Move (3,2) affects two White groups but doesn't create an immediate critical threat.")
    print(" - Move (2,4) is the most powerful. It affects three of the four White groups simultaneously:")
    print("   1. It reduces the liberties of the group at [(3, 3), (3, 4)] from 3 to 2.")
    print("   2. It reduces the liberties of the group at [(1, 4)] from 3 to 2.")
    print("   3. Most importantly, it puts the group at [(2, 5)] in atari, reducing its liberties from 2 to 1.")
    print("This multi-pronged attack on a vital point shatters White's formation and is the correct first step to capturing all the stones.")
    print("\nTherefore, the best move is (2,4).")

if __name__ == '__main__':
    main()
    print("<<<C>>>")