import collections

def get_group_and_liberties(start_stone, all_stones_of_color, opponent_stones, board_size=19):
    """
    Finds a connected group of stones and its liberties using BFS.
    
    Args:
        start_stone (tuple): The (row, col) of the stone to start the search from.
        all_stones_of_color (set): A set of all stones of the same color as start_stone.
        opponent_stones (set): A set of the opponent's stones.
        board_size (int): The size of the board grid.
        
    Returns:
        tuple: A tuple containing:
            - set: The stones belonging to the found group.
            - set: The liberties of that group.
    """
    if start_stone not in all_stones_of_color:
        return set(), set()

    q = collections.deque([start_stone])
    visited = {start_stone}
    group = set()
    liberties = set()

    while q:
        stone = q.popleft()
        group.add(stone)
        
        r, c = stone
        # Check orthogonal neighbors
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            neighbor = (r + dr, c + dc)
            
            # Skip if off-board
            if not (1 <= neighbor[0] <= board_size and 1 <= neighbor[1] <= board_size):
                continue
            
            if neighbor in opponent_stones:
                continue # Neighbor is an opponent stone, not a liberty
            elif neighbor in all_stones_of_color:
                # Neighbor is a friendly stone, add to group if not visited
                if neighbor not in visited:
                    visited.add(neighbor)
                    q.append(neighbor)
            else:
                # Neighbor is an empty point, it's a liberty
                liberties.add(neighbor)
                
    return group, liberties

def analyze_situation(black_stones, white_stones):
    """Analyzes the board to find all white groups and their liberties."""
    
    white_stones_to_check = white_stones.copy()
    groups = []
    
    while white_stones_to_check:
        start_stone = white_stones_to_check.pop()
        group, liberties = get_group_and_liberties(start_stone, white_stones, black_stones)
        
        # Remove all stones of the found group from the set to be checked
        white_stones_to_check -= group
        groups.append({'stones': sorted(list(group)), 'liberties': len(liberties)})
        
    return groups

def main():
    """Main function to solve the Go problem."""
    
    # Initial configuration of stones
    black_stones = {(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)}
    white_stones = {(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)}
    
    # Answer choices (potential moves for Black)
    candidate_moves = {
        "B": (1, 6),
        "C": (2, 4),
        "D": (1, 3),
        "E": (1, 2),
        "F": (3, 2),
        "G": (2, 1),
    }

    print("Initial Analysis:")
    initial_groups = analyze_situation(black_stones, white_stones)
    print(f"Initially, White has {len(initial_groups)} group with {initial_groups[0]['liberties']} liberties.\n")

    best_move = None
    best_move_assessment = ""

    print("--- Analyzing Potential Moves for Black ---")
    for choice, move in candidate_moves.items():
        print(f"\nAnalyzing move {choice}: Black plays at {move}")
        
        # Temporarily add the new black stone
        temp_black_stones = black_stones.union({move})
        
        # Analyze the new board state for White
        new_groups = analyze_situation(temp_black_stones, white_stones)
        
        print(f"  Result: White is split into {len(new_groups)} group(s).")
        for i, group_info in enumerate(new_groups):
            print(f"  - Group {i+1} has {group_info['liberties']} liberties.")
        
        if len(new_groups) > 1:
            assessment = "Excellent. Splitting the group creates multiple weak targets."
            best_move = choice
            best_move_assessment = "This is the vital 'cutting' point. It splits the White stones into two separate, vulnerable groups. White cannot defend both groups at once, so Black can capture them one by one. This move leads to the elimination of all White stones."
        else:
            assessment = "Not effective. The White group remains connected and can likely secure life."
        print(f"  Assessment: {assessment}")

    print("\n--- Conclusion ---")
    if best_move:
        print(f"The best move is choice {best_move}, playing at {candidate_moves[best_move]}.")
        print(best_move_assessment)
    else:
        print("None of the available moves are effective in capturing the white stones.")

if __name__ == "__main__":
    main()