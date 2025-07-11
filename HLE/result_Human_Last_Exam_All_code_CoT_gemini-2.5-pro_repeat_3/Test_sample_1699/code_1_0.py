import collections

# Helper function to find a group of stones and its liberties
def find_group(r, c, board, color):
    """
    Finds a connected group of stones and its liberties.
    
    Args:
        r (int): Starting row of a stone in the group.
        c (int): Starting col of a stone in the group.
        board (dict): The board state.
        color (str): The color of the stones in the group ('B' or 'W').
        
    Returns:
        tuple: A set of stone coordinates in the group and a set of their liberties.
    """
    if (r, c) not in board or board.get((r,c)) != color:
        return set(), set()

    q = collections.deque([(r, c)])
    stones_in_group = set([(r, c)])
    liberties = set()
    visited = set([(r, c)])

    while q:
        row, col = q.popleft()
        
        # Check neighbors (up, down, left, right)
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = row + dr, col + dc
            
            # Skip out-of-bounds coordinates
            if not (1 <= nr <= 19 and 1 <= nc <= 19):
                continue

            neighbor_pos = (nr, nc)
            if neighbor_pos in visited:
                continue
            
            visited.add(neighbor_pos)
            
            if board.get(neighbor_pos) is None: # Empty point
                liberties.add(neighbor_pos)
            elif board.get(neighbor_pos) == color: # Same color stone
                stones_in_group.add(neighbor_pos)
                q.append(neighbor_pos)

    return stones_in_group, liberties

def solve_go_problem():
    """
    Solves the Go problem by analyzing the effect of the move (2, 4).
    """
    black_stones_pos = [(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)]
    white_stones_pos = [(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)]

    # 1. Setup initial board
    board = {}
    for r, c in black_stones_pos:
        board[(r, c)] = 'B'
    for r, c in white_stones_pos:
        board[(r, c)] = 'W'

    # 2. Identify key white groups and their initial liberties
    w_group_1_stones, w_group_1_libs_before = find_group(2, 5, board, 'W')
    w_group_2_stones, w_group_2_libs_before = find_group(3, 4, board, 'W')

    print("Analysis of Black's move at (2, 4):")
    print("="*40)
    
    # 3. Simulate Black's move at (2,4)
    black_move = (2, 4)
    board[black_move] = 'B'

    # 4. Analyze the board state after the move
    _, w_group_1_libs_after = find_group(2, 5, board, 'W')
    _, w_group_2_libs_after = find_group(3, 4, board, 'W')
    
    # 5. Print the impact of the move
    print(f"Impact on White group at {w_group_1_stones}:")
    print(f"Initial liberties: {len(w_group_1_libs_before)}. After B{black_move}, one liberty is taken.")
    print(f"Equation: {len(w_group_1_libs_before)} - 1 = {len(w_group_1_libs_after)}")
    print(f"This group is now in 'atari' (has only {len(w_group_1_libs_after)} liberty).\n")

    print(f"Impact on White group at {w_group_2_stones}:")
    print(f"Initial liberties: {len(w_group_2_libs_before)}. After B{black_move}, one liberty is taken.")
    print(f"Equation: {len(w_group_2_libs_before)} - 1 = {len(w_group_2_libs_after)}")
    print(f"This group is now severely weakened with only {len(w_group_2_libs_after)} liberties.\n")
    
    print("Conclusion:")
    print(f"The move at (2, 4) is the critical move (tesuji). It attacks two groups simultaneously,")
    print("putting one in immediate danger (atari) and restricting the other. This creates a")
    print("situation where White cannot defend all of its stones, allowing Black to eventually")
    print("capture all of them. Therefore, the correct first move is (2, 4).")


solve_go_problem()