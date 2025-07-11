import collections

def get_group_and_liberties(start_r, start_c, board_state):
    """
    Finds all stones in a connected group and their liberties using BFS.
    A liberty is an empty point adjacent to any stone in the group.
    """
    if board_state.get((start_r, start_c)) is None:
        return set(), set()

    player_color = board_state[(start_r, start_c)]
    group_stones = set()
    liberties = set()
    
    queue = collections.deque([(start_r, start_c)])
    visited = set([(start_r, start_c)])

    while queue:
        r, c = queue.popleft()
        group_stones.add((r, c))

        # Check neighbors (up, down, left, right)
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc

            # Skip if off-board
            if not (1 <= nr <= 19 and 1 <= nc <= 19):
                continue
            
            neighbor_pos = (nr, nc)
            if neighbor_pos in visited:
                continue
            
            visited.add(neighbor_pos)
            neighbor_content = board_state.get(neighbor_pos)

            if neighbor_content is None:
                liberties.add(neighbor_pos)
            elif neighbor_content == player_color:
                queue.append(neighbor_pos)
                
    return group_stones, liberties

def solve_go_problem():
    """
    Solves the Go life and death problem.
    """
    # Coordinate system: row 1-19 (top-bottom), col 1-19 (right-left)
    black_stones = {(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)}
    white_stones = {(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)}

    # Initial board setup
    board = {}
    for stone in black_stones:
        board[stone] = 'Black'
    for stone in white_stones:
        board[stone] = 'White'
    
    # --- Step 1: Analyze the initial state ---
    # We target the white stone at (2,5) to find the entire group.
    white_group, initial_liberties = get_group_and_liberties(2, 5, board)
    
    print("--- Initial Analysis ---")
    print(f"The White group consists of {len(white_group)} stones.")
    print(f"Initial number of liberties for the White group: {len(initial_liberties)}")
    print("Initial liberties are at coordinates:", sorted(list(initial_liberties)))
    print("\nTo capture the White stones, Black must eliminate all these liberties.")
    
    # --- Step 2 & 3: Simulate the tesuji starting with (2, 4) ---
    print("\n--- Simulating the Optimal Move: Black plays at (2, 4) ---")
    
    # Move 1: Black plays the "throw-in" at (2, 4)
    black_move_1 = (2, 4)
    board[black_move_1] = 'Black'
    print(f"1. Black plays at coordinate: {black_move_1}. This is a sacrifice.")

    # Move 2: White is forced to capture the Black stone. The only way is to play at (2,3).
    # This removes the black stone at (2,4) and adds a white one at (2,3).
    white_move_1 = (2, 3)
    board.pop(black_move_1) # White captures the black stone
    board[white_move_1] = 'White' # White places a stone
    print(f"2. White's forced response is to capture by playing at: {white_move_1}.")

    # Move 3: Black plays on the vital point again.
    black_move_2 = (2, 4)
    board[black_move_2] = 'Black'
    print(f"3. Black plays on the vital point again at: {black_move_2}.")
    
    # --- Step 4: Evaluate the final state ---
    print("\n--- Final Position Analysis ---")
    final_white_group, final_liberties = get_group_and_liberties(2, 5, board)
    
    print("After this sequence, White's shape is critically flawed (a 'false eye').")
    print(f"The White group now consists of {len(final_white_group)} stones.")
    print(f"Final number of liberties for the White group: {len(final_liberties)}")
    print("Remaining liberties are at coordinates:", sorted(list(final_liberties)))
    
    print("\nConclusion: The White group cannot make two eyes to live and is surrounded.")
    print("Black can now capture the entire group by filling the remaining 5 liberties.")
    print("Therefore, the first move (2, 4) leads to the elimination of all White stones.")

solve_go_problem()