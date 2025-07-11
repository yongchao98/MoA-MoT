def solve_go_problem():
    """
    Analyzes the Go board to find the optimal move for Black.
    """
    board_str = [
        "EWBEEEBWE",
        "EWBEEWBWE",
        "EWWBEWBWE",
        "EEWBEEBWE",
        "EEWWBBBWE",
        "EEEEWWWWE",
        "EEEEEEEEE",
        "EEEEEEEEE",
    ]

    board = [list(row) for row in board_str]
    rows, cols = len(board), len(board[0])

    def get_neighbors(r, c):
        """Returns valid neighbors for a given coordinate."""
        neighbors = []
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                neighbors.append((nr, nc))
        return neighbors

    def find_groups(player):
        """Finds all groups and their liberties for a given player."""
        groups = []
        visited = set()
        for r in range(rows):
            for c in range(cols):
                if board[r][c] == player and (r, c) not in visited:
                    group_stones = set()
                    group_liberties = set()
                    q = [(r, c)]
                    visited.add((r, c))
                    
                    head = 0
                    while head < len(q):
                        curr_r, curr_c = q[head]
                        head += 1
                        group_stones.add((curr_r, curr_c))
                        for nr, nc in get_neighbors(curr_r, curr_c):
                            if board[nr][nc] == player and (nr, nc) not in visited:
                                visited.add((nr, nc))
                                q.append((nr, nc))
                            elif board[nr][nc] == 'E':
                                group_liberties.add((nr, nc))
                    groups.append({'stones': group_stones, 'liberties': group_liberties})
        return groups

    # 1. Find all white and black groups
    white_groups = find_groups('W')
    black_groups = find_groups('B')

    # 2. Identify the target group (2 white stones) and the black groups
    target_white_group = next((g for g in white_groups if len(g['stones']) == 2), None)
    
    # Sort black groups by number of liberties to find the weak and strong ones
    black_groups.sort(key=lambda g: len(g['liberties']))
    weak_black_group = black_groups[0]
    strong_black_group = black_groups[1]
    
    # 3. Find the connecting point between the two black groups
    # This is a shared liberty, the key to survival.
    connecting_points = weak_black_group['liberties'].intersection(strong_black_group['liberties'])

    solutions = []
    # 4. Analyze the connecting move
    for move in connecting_points:
        # Check if this move wins the capturing race
        
        # Calculate liberties after the move
        # Black's new liberties = all libs of both groups, minus the move itself, plus any new ones.
        new_black_libs = weak_black_group['liberties'].union(strong_black_group['liberties'])
        new_black_libs.remove(move)
        for nr, nc in get_neighbors(move[0], move[1]):
            if board[nr][nc] == 'E':
                 new_black_libs.add((nr, nc))

        # White's liberties are unaffected by this move as it's not adjacent
        white_libs_count = len(target_white_group['liberties'])
        
        # In a race (semeai), after Black moves, it is White's turn.
        # White needs to fill black_libs_count liberties. Black needs to fill white_libs_count.
        # Black wins if they need fewer moves.
        if white_libs_count < len(new_black_libs):
            # The move is a solution. Convert 0-indexed to 1-indexed for output.
            solutions.append((move[0] + 1, move[1] + 1))
    
    # Sort lexicographically and print
    solutions.sort()
    # Format the output string like (r, c), (r, c)
    print(", ".join(f"({r}, {c})" for r, c in solutions))

solve_go_problem()