from collections import deque

def get_knight_moves(pos):
    row, col = pos
    moves = [
        (row+2, col+1), (row+2, col-1),
        (row-2, col+1), (row-2, col-1),
        (row+1, col+2), (row+1, col-2),
        (row-1, col+2), (row-1, col-2)
    ]
    return [(r, c) for r, c in moves if 0 <= r < 3 and 0 <= c < 4]

def pos_to_notation(pos):
    return f"{chr(65+pos[1])}{3-pos[0]}"

def evaluate_position(whites, blacks, target_w, target_b):
    # Calculate Manhattan distance to targets
    total_distance = 0
    for w in whites:
        min_dist = min(abs(w[0]-t[0]) + abs(w[1]-t[1]) for t in target_w)
        total_distance += min_dist
    for b in blacks:
        min_dist = min(abs(b[0]-t[0]) + abs(b[1]-t[1]) for t in target_b)
        total_distance += min_dist
    return total_distance

def find_optimal_solution():
    initial_board = [
        ['.', 'w', 'B', ' '],  # row 3
        ['B', ' ', ' ', '.'],  # row 2
        [' ', '.', 'w', ' ']   # row 1
    ]
    
    target_w = {(0, 2), (1, 0)}  # Target positions for white (original black positions)
    target_b = {(0, 1), (2, 2)}  # Target positions for black (original white positions)
    
    queue = deque([(initial_board, [], True, 0)])  # board, moves, is_white_turn, depth
    visited = {str(initial_board): 0}  # state: minimum depth seen
    max_depth = 8  # Reasonable limit for optimal solution
    
    while queue:
        current_board, moves, is_white_turn, depth = queue.popleft()
        
        if depth >= max_depth:
            continue
        
        # Get current positions
        whites = []
        blacks = []
        empty = []
        for i in range(3):
            for j in range(4):
                if current_board[i][j] == 'w':
                    whites.append((i, j))
                elif current_board[i][j] == 'B':
                    blacks.append((i, j))
                elif current_board[i][j] == '.':
                    empty.append((i, j))
        
        # Check if we reached target state
        if set(whites) == target_w and set(blacks) == target_b:
            return moves
        
        # Generate moves for current player
        current_pieces = whites if is_white_turn else blacks
        piece_symbol = 'w' if is_white_turn else 'B'
        
        # Sort moves by evaluation score
        possible_moves = []
        for piece_pos in current_pieces:
            for new_pos in get_knight_moves(piece_pos):
                if new_pos in empty:
                    new_board = [row[:] for row in current_board]
                    new_board[piece_pos[0]][piece_pos[1]] = '.'
                    new_board[new_pos[0]][new_pos[1]] = piece_symbol
                    
                    # Calculate new positions for evaluation
                    new_whites = [(i, j) for i in range(3) for j in range(4) if new_board[i][j] == 'w']
                    new_blacks = [(i, j) for i in range(3) for j in range(4) if new_board[i][j] == 'B']
                    
                    score = evaluate_position(new_whites, new_blacks, target_w, target_b)
                    possible_moves.append((score, piece_pos, new_pos, new_board))
        
        # Sort moves by evaluation score
        possible_moves.sort(key=lambda x: x[0])
        
        # Try the most promising moves first
        for _, piece_pos, new_pos, new_board in possible_moves:
            board_str = str(new_board)
            if board_str not in visited or visited[board_str] > depth + 1:
                visited[board_str] = depth + 1
                new_moves = moves + [f"{piece_symbol},{pos_to_notation(piece_pos)},{pos_to_notation(new_pos)}"]
                queue.append((new_board, new_moves, not is_white_turn, depth + 1))
    
    return None

solution = find_optimal_solution()
if solution:
    print(solution)
else:
    print("No")