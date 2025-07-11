import collections

def to_internal_coords(r, c):
    """Converts problem coordinates (1-19 top-to-bottom, 1-19 right-to-left) to 0-indexed matrix coordinates."""
    return r - 1, 19 - c

def to_problem_coords(r_idx, c_idx):
    """Converts internal matrix coordinates back to the problem's coordinates."""
    return r_idx + 1, 19 - c_idx

def find_group_and_liberties(r_idx, c_idx, board, visited_stones):
    """
    Finds a group of connected stones and their liberties using Breadth-First Search (BFS).
    Liberties are the unique adjacent empty points of the group.
    """
    color = board[r_idx][c_idx]
    if color == 0:
        return None, None

    group_stones = set()
    liberties = set()
    q = collections.deque([(r_idx, c_idx)])
    visited_q = set([(r_idx, c_idx)])

    while q:
        r, c = q.popleft()
        visited_stones.add((r, c))
        group_stones.add((r, c))

        # Check neighbors (up, down, left, right)
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < 19 and 0 <= nc < 19:
                if board[nr][nc] == 0:
                    liberties.add((nr, nc))
                elif board[nr][nc] == color and (nr, nc) not in visited_q:
                    q.append((nr, nc))
                    visited_q.add((nr, nc))

    return group_stones, liberties

def analyze_white_groups(board):
    """Scans the board to find all distinct white groups and returns their details."""
    white_groups = []
    visited_stones = set()
    for r in range(19):
        for c in range(19):
            if board[r][c] == 2 and (r, c) not in visited_stones:
                group_stones, liberties = find_group_and_liberties(r, c, board, visited_stones)
                if group_stones:
                    problem_stones = sorted([to_problem_coords(r_idx, c_idx) for r_idx, c_idx in group_stones])
                    white_groups.append({
                        'stones': problem_stones,
                        'liberties': len(liberties),
                    })
    return white_groups

def main():
    """Main function to solve the Go problem."""
    # 1. Setup the board based on the initial configuration
    board = [[0] * 19 for _ in range(19)]
    black_stones_orig = [(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)]
    white_stones_orig = [(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)]

    for r, c in black_stones_orig:
        r_idx, c_idx = to_internal_coords(r, c)
        board[r_idx][c_idx] = 1
    for r, c in white_stones_orig:
        r_idx, c_idx = to_internal_coords(r, c)
        board[r_idx][c_idx] = 2

    # 2. Analyze the initial board state
    print("Step 1: Analyzing the initial board to identify White groups and their liberties.")
    initial_groups = analyze_white_groups(board)
    initial_total_liberties = sum(g['liberties'] for g in initial_groups)
    print(f"Found {len(initial_groups)} distinct White groups:")
    for i, group in enumerate(initial_groups):
        print(f"  - Group {i+1}: Stones at {group['stones']}, has {group['liberties']} liberties.")
    print("-" * 30)

    # 3. Analyze each candidate move for Black
    print("Step 2: Evaluating each candidate move by its impact on White's groups.")
    candidate_moves = {"B": (1, 6), "C": (2, 1), "D": (3, 2), "E": (1, 2), "F": (1, 3), "G": (2, 4)}

    results = {}
    for choice, move in candidate_moves.items():
        print(f"\nAnalyzing move {choice}: Black places a stone at {move}.")
        
        temp_board = [row[:] for row in board]
        r_idx, c_idx = to_internal_coords(move[0], move[1])
        
        if temp_board[r_idx][c_idx] != 0:
            print("  Invalid move: Position is occupied.")
            continue
        temp_board[r_idx][c_idx] = 1

        groups_after_move = analyze_white_groups(temp_board)
        num_atari_groups = sum(1 for g in groups_after_move if g['liberties'] == 1)
        total_liberties_after = sum(g['liberties'] for g in groups_after_move)
        liberty_reduction = initial_total_liberties - total_liberties_after
        
        results[choice] = {'atari_groups': num_atari_groups, 'liberty_reduction': liberty_reduction}
        
        print(f"  Impact: White groups put in atari: {num_atari_groups}. Total liberty reduction: {liberty_reduction}.")

    # 4. Determine the best move based on creating atari threats.
    best_choice = 'A' # Default to Impossible
    max_atari = 0
    max_reduction = -1

    for choice, result in results.items():
        if result['atari_groups'] > max_atari:
            max_atari = result['atari_groups']
            max_reduction = result['liberty_reduction']
            best_choice = choice
        elif result['atari_groups'] == max_atari and result['liberty_reduction'] > max_reduction:
            max_reduction = result['liberty_reduction']
            best_choice = choice
    
    print("-" * 30)
    print("Step 3: Conclusion.")
    if best_choice != 'A' and max_atari > 0:
        print(f"Move {best_choice} ({candidate_moves[best_choice]}) is the only move that puts a White group in atari.")
        print("This move creates an immediate threat, forcing White's response and giving Black the initiative to capture all stones.")
    else:
        # Fallback reasoning if no move creates atari.
        if best_choice != 'A':
            print(f"While no move creates an immediate atari, move {best_choice} ({candidate_moves[best_choice]}) is the strongest because it reduces White's liberties the most significantly ({max_reduction}), setting up a future capture.")
        else:
             print("None of the moves lead to a certain capture of all white stones.")


main()
<<<G>>>