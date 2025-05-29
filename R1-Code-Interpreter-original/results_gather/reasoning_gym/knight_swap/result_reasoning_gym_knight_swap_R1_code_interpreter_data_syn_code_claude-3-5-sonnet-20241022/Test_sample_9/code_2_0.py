from collections import deque

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(x, y) for x, y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return tuple(tuple(row) for row in board)

def is_goal_state(board):
    # White knights should be at B3 and C1 (original black positions)
    # Black knights should be at B1 and D2 (original white positions)
    return (board[0][1] == 'w' and board[2][2] == 'w' and 
            board[2][1] == 'B' and board[1][3] == 'B')

def make_move(board, start, end):
    new_board = [list(row) for row in board]
    piece = new_board[start[1]][start[0]]
    new_board[start[1]][start[0]] = '.'
    new_board[end[1]][end[0]] = piece
    return [list(row) for row in new_board]

def solve_puzzle():
    # Initial board state (0-based indexing)
    initial_board = [
        ['.', 'B', '.', '.'],  # row 3
        ['.', '.', '.', 'w'],  # row 2
        ['.', 'w', 'B', '.']   # row 1
    ]
    
    # Starting positions
    white_initial = [(1,2), (3,1)]  # B1, D2
    black_initial = [(1,0), (2,2)]  # B3, C1
    
    queue = deque([(initial_board, [], True, white_initial, black_initial)])
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn, white_pos, black_pos = queue.popleft()
        
        if len(moves) > 20:  # Prevent infinite loops
            continue
            
        if is_goal_state(board):
            return moves
        
        current_positions = white_pos if is_white_turn else black_pos
        piece = 'w' if is_white_turn else 'B'
        
        for start_pos in current_positions:
            for end_pos in get_knight_moves(start_pos):
                if board[end_pos[1]][end_pos[0]] == '.':
                    new_board = make_move(board, start_pos, end_pos)
                    board_state = board_to_string(new_board)
                    
                    if board_state not in visited:
                        visited.add(board_state)
                        
                        # Update positions
                        new_white_pos = [(end_pos[0], end_pos[1]) if pos == start_pos else pos for pos in white_pos] if is_white_turn else white_pos
                        new_black_pos = [(end_pos[0], end_pos[1]) if pos == start_pos else pos for pos in black_pos] if not is_white_turn else black_pos
                        
                        # Record move
                        move = f"{piece},{chr(65+start_pos[0])}{start_pos[1]+1},{chr(65+end_pos[0])}{end_pos[1]+1}"
                        queue.append((new_board, moves + [move], not is_white_turn, new_white_pos, new_black_pos))
    
    return None

solution = solve_puzzle()
if solution:
    print(solution)
else:
    print("No")