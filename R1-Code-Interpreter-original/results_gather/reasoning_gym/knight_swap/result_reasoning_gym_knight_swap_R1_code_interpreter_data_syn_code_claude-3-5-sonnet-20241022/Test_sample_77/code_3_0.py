from collections import deque

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(str(cell) for cell in row) for row in board)

def get_positions(board, piece):
    return [(i,j) for i in range(4) for j in range(3) 
            if board[j][i] == piece]

def make_move(board, start, end):
    new_board = [row[:] for row in board]
    piece = new_board[start[1]][start[0]]
    new_board[start[1]][start[0]] = '.'
    new_board[end[1]][end[0]] = piece
    return new_board

def is_goal_state(board):
    # Check if white knights are at original black positions (A3, C2)
    # and black knights are at original white positions (B1, D2)
    white_positions = set(get_positions(board, 'w'))
    black_positions = set(get_positions(board, 'B'))
    
    target_white = {(0,0), (2,1)}  # A3, C2
    target_black = {(1,2), (3,1)}  # B1, D2
    
    return white_positions == target_white and black_positions == target_black

def solve_knights():
    initial_board = [
        ['B', '.', '.', ' '],  # row 3 (index 0)
        [' ', ' ', 'B', 'w'],  # row 2 (index 1)
        ['.', 'w', '.', '.']   # row 1 (index 2)
    ]
    
    queue = deque([(initial_board, [], True)])
    visited = {(board_to_string(initial_board), True)}
    
    while queue:
        current_board, moves, is_white_turn = queue.popleft()
        
        if is_goal_state(current_board):
            return moves
        
        # Get current player's pieces
        current_pieces = get_positions(current_board, 'w' if is_white_turn else 'B')
        
        for start in current_pieces:
            for end in get_knight_moves(start):
                if current_board[end[1]][end[0]] == '.':
                    new_board = make_move(current_board, start, end)
                    state = (board_to_string(new_board), not is_white_turn)
                    
                    if state not in visited:
                        visited.add(state)
                        # Convert coordinates to chess notation
                        move = f"{'w' if is_white_turn else 'B'},{chr(65+start[0])}{3-start[1]},{chr(65+end[0])}{3-end[1]}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knights()
if solution:
    print(solution)
else:
    print("No")