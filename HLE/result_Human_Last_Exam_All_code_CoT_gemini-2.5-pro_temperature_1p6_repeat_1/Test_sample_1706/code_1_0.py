import collections

# Board configuration from the problem description
BLACK_STONES = {(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)}
WHITE_STONES = {(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)}
CANDIDATE_MOVES = {
    "B": (1, 6),
    "C": (2, 1),
    "D": (3, 2),
    "E": (1, 2),
    "F": (1, 3),
    "G": (2, 4)
}
BOARD_SIZE = 19

def get_neighbors(r, c):
    """Returns a list of valid neighbor coordinates for a given point."""
    neighbors = []
    if r > 1: neighbors.append((r - 1, c))
    if r < BOARD_SIZE: neighbors.append((r + 1, c))
    if c > 1: neighbors.append((r, c - 1))
    if c < BOARD_SIZE: neighbors.append((r, c + 1))
    return neighbors

def find_group_and_liberties(start_r, start_c, color_stones, all_stones):
    """
    Finds a connected group of stones and its liberties using Breadth-First Search.
    
    Args:
        start_r, start_c: The coordinates of a stone to start the search from.
        color_stones: A set of all stones of the same color as the starting stone.
        all_stones: A set of all stones on the board (black and white).

    Returns:
        A tuple containing:
        - A set of (r, c) tuples representing the stones in the group.
        - A set of (r, c) tuples representing the liberties of the group.
    """
    if (start_r, start_c) not in color_stones:
        return set(), set()

    q = collections.deque([(start_r, start_c)])
    group = set(q)
    liberties = set()
    
    while q:
        r, c = q.popleft()
        for neighbor in get_neighbors(r, c):
            if neighbor in group or neighbor in liberties:
                continue
            
            if neighbor in color_stones:
                group.add(neighbor)
                q.append(neighbor)
            elif neighbor not in all_stones:
                liberties.add(neighbor)
    
    return group, liberties

def analyze_situation(black_stones, white_stones, candidate_moves):
    """
    Analyzes the board to find the best move for Black from a list of candidates.
    """
    all_stones = black_stones.union(white_stones)
    
    # 1. Find all initial white groups and their liberties
    initial_white_groups = []
    stones_already_in_group = set()
    for r, c in white_stones:
        if (r, c) not in stones_already_in_group:
            group, liberties = find_group_and_liberties(r, c, white_stones, all_stones)
            initial_white_groups.append({'stones': group, 'liberties': liberties})
            stones_already_in_group.update(group)
            
    best_move = None
    best_score = -1
    best_move_calc_str = ""

    # 2. Evaluate each candidate move
    for label, move in candidate_moves.items():
        if move in all_stones:
            continue
        
        atari_created = 0
        libs_removed_count = 0
        
        # Calculate score based on impact on white groups
        for group_info in initial_white_groups:
            if move in group_info['liberties']:
                libs_removed_count += 1
                # If the move takes the last liberty but one, it's an atari
                if len(group_info['liberties']) == 2:
                    atari_created += 1

        # Scoring: 10 points for creating an atari, 1 point for each group affected.
        score = (atari_created * 10) + libs_removed_count
        
        if score > best_score:
            best_score = score
            best_move = move
            best_move_calc_str = f"atari_groups_created * 10 + liberties_removed_from = score\n{atari_created} * 10 + {libs_removed_count} = {score}"

    # 3. Print the result
    print("Analysis complete.")
    print(f"The best move is to play at {best_move}.")
    print("This move is the most effective because it attacks multiple white groups simultaneously.")
    print("\nScoring calculation for the best move:")
    print(best_move_calc_str)
    
# Run the analysis
analyze_situation(BLACK_STONES, WHITE_STONES, CANDIDATE_MOVES)
