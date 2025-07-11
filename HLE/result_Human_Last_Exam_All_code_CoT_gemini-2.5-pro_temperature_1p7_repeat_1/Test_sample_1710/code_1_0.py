def solve_go_problem():
    """
    Solves the Go problem by finding the vital move for Black.

    The logic is to:
    1. Parse the board representation.
    2. Define a function to find connected groups of stones and their liberties.
    3. Identify all black groups on the board.
    4. Find the specific black group that is in 'atari' (has only 2 liberties).
       The liberties of this group are the only candidate moves for survival.
    5. From the candidate moves, determine which one is strategically superior.
       The best move is the one that connects the weak group to another friendly group,
       as this strengthens the overall position and enables the capture described in the prompt.
    6. Format and print the final answer.
    """
    board_str = "EWBEEEBWE EWBEEWBWE EWWBEWBWE EEWBEEBWE EEWWBBBWE EEEEWWWWE EEEEEEEEE EEEEEEEEE"
    board = [list(row) for row in board_str.split()]
    height = len(board)
    width = len(board[0])

    def get_neighbors(r, c):
        neighbors = []
        if r > 0: neighbors.append((r - 1, c))
        if r < height - 1: neighbors.append((r + 1, c))
        if c > 0: neighbors.append((r, c - 1))
        if c < width - 1: neighbors.append((r, c + 1))
        return neighbors

    def find_group(start_r, start_c, player, global_visited):
        if (start_r, start_c) in global_visited or board[start_r][start_c] != player:
            return None, None

        q = [(start_r, start_c)]
        group_visited = set(q)
        group_stones = set(q)
        liberties = set()

        head = 0
        while head < len(q):
            r, c = q[head]
            head += 1
            global_visited.add((r, c))

            for nr, nc in get_neighbors(r, c):
                if (nr, nc) not in group_visited:
                    neighbor_stone = board[nr][nc]
                    if neighbor_stone == player:
                        group_visited.add((nr, nc))
                        q.append((nr, nc))
                        group_stones.add((nr, nc))
                    elif neighbor_stone == 'E':
                        liberties.add((nr, nc))
        return group_stones, liberties

    # Find all black groups and identify the one to be saved
    all_black_groups = []
    visited_stones = set()
    for r in range(height):
        for c in range(width):
            if board[r][c] == 'B' and (r, c) not in visited_stones:
                stones, libs = find_group(r, c, 'B', visited_stones)
                if stones:
                    all_black_groups.append({'stones': stones, 'liberties': libs})

    survival_candidate_moves = []
    weak_group_stones = set()
    for group in all_black_groups:
        if len(group['liberties']) == 2:
            survival_candidate_moves = list(group['liberties'])
            weak_group_stones = group['stones']
            break

    # Find the move that connects the weak group to another black group
    best_move = None
    for move_r, move_c in survival_candidate_moves:
        is_connecting_move = False
        for nr, nc in get_neighbors(move_r, move_c):
            # Check if this neighbor is a black stone NOT in our weak group
            if board[nr][nc] == 'B' and (nr, nc) not in weak_group_stones:
                is_connecting_move = True
                break
        if is_connecting_move:
            best_move = (move_r + 1, move_c + 1) # Convert to 1-based index
            break
            
    # Fallback in case the connecting move logic doesn't find a unique solution,
    # though it does in this case. Based on analysis, (2,4) is the only correct answer.
    if best_move:
        print(f"({best_move[0]}, {best_move[1]})")
    else:
        # Sort lexicographically and print if multiple moves were found without a clear best one.
        # This part of the code is unlikely to be reached for this specific problem.
        results = sorted([(r + 1, c + 1) for r, c in survival_candidate_moves])
        print(', '.join([f"({r}, {c})" for r, c in results]))

solve_go_problem()