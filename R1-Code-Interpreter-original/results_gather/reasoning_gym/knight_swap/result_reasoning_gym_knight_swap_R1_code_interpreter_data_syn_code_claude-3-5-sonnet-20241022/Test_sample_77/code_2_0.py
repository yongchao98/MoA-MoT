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
    # Original black positions: A3, C2
    # Original white positions: B1, D2
    white_pos = set((x,y) for x in range(4) for y in range(3) if board[y][x] == 'w')
    black_pos = set((x,y) for x in range(4) for y in range(3) if board[y][x] == 'B')
    return (0,0) in white_pos and (2,1) in white_pos and (1,2) in black_pos and (3,1) in black_pos

def solve_knights():
    # Initialize board: 0-based indexing
    initial_board = [
        ['B', '.', '.', ' '],  # row 3
        [' ', ' ', 'B', 'w'],  # row 2
        ['.', 'w', '.', '.']   # row 1
    ]
    
    queue = deque([(initial_board, [], True)])  # (board, moves_list, is_white_turn)
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, is_white_turn = queue.popleft()
        
        if is_goal_state(current_board):
            return moves
            
        current_pieces = get_positions(current_board, 'w' if is_white_turn else 'B')
        
        for start in current_pieces:
            for end in get_knight_moves(start):
                if current_board[end[1]][end[0]] == '.':
                    new_board = make_move(current_board, start, end)
                    board_str = board_to_string(new_board)
                    
                    if board_str not in visited:
                        visited.add(board_str)
                        # Convert to chess notation (A1 is bottom-left)
                        move = f"{'w' if is_white_turn else 'B'},{chr(65+start[0])}{3-start[1]},{chr(65+end[0])}{3-end[1]}"
                        queue.append((new_board, moves + [move], not is_white_turn))

    return None

solution = solve_knights()
if solution:
    print(solution)
else:
    print("No")