from collections import deque

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def pos_to_coord(pos):
    return (ord(pos[0]) - ord('A'), int(pos[1]) - 1)

def coord_to_pos(coord):
    return chr(coord[0] + ord('A')) + str(coord[1] + 1)

def is_target_state(board):
    # White knights should be at B3,C3 (originally black positions)
    # Black knights should be at A2,D2 (originally white positions)
    return (board[0][1] == 'w' and board[0][2] == 'w' and
            board[1][0] == 'B' and board[1][3] == 'B')

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def solve_knight_swap():
    # Initial state: Black knights at B3,C3; White knights at A2,D2
    initial_board = [
        ['.', 'B', 'B', '.'],  # row 3
        ['w', '.', '.', 'w'],  # row 2
        ['.', '.', '.', '.']   # row 1
    ]
    
    queue = deque([(initial_board, [], True)])  # (board, moves, is_white_turn)
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, is_white_turn = queue.popleft()
        current_piece = 'w' if is_white_turn else 'B'
        
        # Find all knights of current color
        knights = []
        for y in range(3):
            for x in range(4):
                if current_board[y][x] == current_piece:
                    knights.append((x, y))
        
        # Try each knight's possible moves
        for kx, ky in knights:
            for new_x, new_y in get_knight_moves((kx, ky)):
                if current_board[new_y][new_x] == '.':  # Only move to empty squares
                    # Create new board state
                    new_board = [row[:] for row in current_board]
                    new_board[ky][kx] = '.'
                    new_board[new_y][new_x] = current_piece
                    
                    board_state = board_to_string(new_board)
                    if board_state not in visited:
                        visited.add(board_state)
                        from_pos = coord_to_pos((kx, ky))
                        to_pos = coord_to_pos((new_x, new_y))
                        new_moves = moves + [f"{current_piece},{from_pos},{to_pos}"]
                        
                        if is_target_state(new_board):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return "No"

solution = solve_knight_swap()
print(solution)