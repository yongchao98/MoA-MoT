from collections import deque

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_state(board):
    # Convert board to immutable state including piece positions
    w_pos = tuple(sorted((i,j) for i in range(4) for j in range(3) if board[j][i] == 'w'))
    b_pos = tuple(sorted((i,j) for i in range(4) for j in range(3) if board[j][i] == 'B'))
    return (w_pos, b_pos)

def make_move(board, start, end):
    new_board = [[cell for cell in row] for row in board]
    piece = new_board[start[1]][start[0]]
    new_board[start[1]][start[0]] = '.'
    new_board[end[1]][end[0]] = piece
    return new_board

def is_goal_state(white_pos, black_pos):
    target_white = {(0,0), (2,1)}  # A3, C2
    target_black = {(1,2), (3,1)}  # B1, D2
    return set(white_pos) == target_white and set(black_pos) == target_black

def solve_knights():
    initial_board = [
        ['B', '.', '.', ' '],  # row 3
        [' ', ' ', 'B', 'w'],  # row 2
        ['.', 'w', '.', '.']   # row 1
    ]
    
    # Initial positions
    start_state = board_to_state(initial_board)
    queue = deque([(initial_board, [], True, start_state)])
    visited = {(start_state, True)}
    
    while queue:
        current_board, moves, is_white_turn, current_state = queue.popleft()
        
        if len(moves) > 8:  # Limit depth as optimal solution should be shorter
            continue
            
        white_pos, black_pos = current_state
        if is_goal_state(white_pos, black_pos):
            return moves
        
        # Get current player's pieces and their positions
        current_pieces = [(i,j) for i,j in (white_pos if is_white_turn else black_pos)]
        
        for start in current_pieces:
            for end in get_knight_moves(start):
                if current_board[end[1]][end[0]] == '.':
                    new_board = make_move(current_board, start, end)
                    new_state = board_to_state(new_board)
                    state_key = (new_state, not is_white_turn)
                    
                    if state_key not in visited:
                        visited.add(state_key)
                        move = f"{'w' if is_white_turn else 'B'},{chr(65+start[0])}{3-start[1]},{chr(65+end[0])}{3-end[1]}"
                        queue.append((new_board, moves + [move], not is_white_turn, new_state))
    
    return None

solution = solve_knights()
if solution:
    print(solution)
else:
    print("No")