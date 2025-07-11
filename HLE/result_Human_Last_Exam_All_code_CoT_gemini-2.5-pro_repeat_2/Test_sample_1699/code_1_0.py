import collections

def get_neighbors(r, c):
    """Returns a list of valid neighbor coordinates."""
    neighbors = []
    if r > 1: neighbors.append((r - 1, c))
    if r < 19: neighbors.append((r + 1, c))
    if c > 1: neighbors.append((r, c - 1))
    if c < 19: neighbors.append((r, c + 1))
    return neighbors

def find_group(start_node, all_stones):
    """Finds all connected stones of the same color (a group)."""
    if start_node not in all_stones:
        return set(), set()
    
    q = collections.deque([start_node])
    visited = {start_node}
    group = {start_node}
    
    while q:
        curr_r, curr_c = q.popleft()
        for neighbor in get_neighbors(curr_r, curr_c):
            if neighbor in all_stones and neighbor not in visited:
                visited.add(neighbor)
                group.add(neighbor)
                q.append(neighbor)
    return group

def get_liberties(group, black_stones, white_stones):
    """Calculates the set of liberties for a given group."""
    liberties = set()
    for r, c in group:
        for neighbor in get_neighbors(r, c):
            if neighbor not in black_stones and neighbor not in white_stones:
                liberties.add(neighbor)
    return liberties

def solve_go_problem():
    """
    Analyzes the Go board state to find the best move for Black.
    """
    black_stones = {(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)}
    white_stones = {(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)}
    
    # All white stones form a single group
    white_start_stone = next(iter(white_stones))
    white_group = find_group(white_start_stone, white_stones)
    
    initial_liberties = get_liberties(white_group, black_stones, white_stones)
    print(f"The White group initially has {len(initial_liberties)} liberties: {sorted(list(initial_liberties))}\n")

    # Answer choices for Black's move
    candidate_moves = {
        "C": (2, 4),
        "D": (1, 3),
        "E": (1, 2),
        "F": (3, 2),
        "G": (2, 1),
    }

    best_black_move = None
    min_max_libs_for_white = float('inf')

    print("Analyzing Black's candidate moves...")
    for choice, b_move in candidate_moves.items():
        # Step 1: Black makes a move
        temp_black_stones = black_stones | {b_move}
        
        # The potential responses for White are the liberties of its original group
        white_responses = get_liberties(white_group, black_stones, white_stones)
        if b_move in white_responses:
            white_responses.remove(b_move)
        
        max_libs_for_white = -1
        best_white_response = None
        
        # Step 2: Find White's best response (the one that maximizes its liberties)
        for w_move in white_responses:
            temp_white_stones = white_stones | {w_move}
            
            # Find the new white group after White's move
            new_white_group = find_group(w_move, temp_white_stones)
            
            # Calculate liberties of this new group
            libs = get_liberties(new_white_group, temp_black_stones, temp_white_stones)
            
            if len(libs) > max_libs_for_white:
                max_libs_for_white = len(libs)
                best_white_response = w_move

        print(f"Move {choice} {b_move}:")
        print(f"  If Black plays at {b_move}, White's best response is {best_white_response},")
        print(f"  which leaves the White group with {max_libs_for_white} liberties.\n")

        if max_libs_for_white < min_max_libs_for_white:
            min_max_libs_for_white = max_libs_for_white
            best_black_move = (choice, b_move)
            
    print("--- Conclusion ---")
    print("To capture the White group, Black must choose the move that is most severe.")
    print("This means picking the move that leaves White with the minimum number of liberties, even after White makes its best defensive reply.")
    print(f"The move that minimizes White's maximum achievable liberties is {best_black_move[1]}.")
    print(f"This move is choice {best_black_move[0]}.")

# Run the analysis
solve_go_problem()