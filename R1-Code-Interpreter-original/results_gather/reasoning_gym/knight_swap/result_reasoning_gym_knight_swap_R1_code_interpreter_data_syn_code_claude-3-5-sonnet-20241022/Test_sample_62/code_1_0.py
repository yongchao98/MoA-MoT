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

def pos_to_coord(pos):
    return chr(pos[1] + ord('A')) + str(3 - pos[0])  # Corrected coordinate system

def coord_to_pos(coord):
    return (3 - int(coord[1]), ord(coord[0]) - ord('A'))

def get_board_state(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # White knights should be at B1, C1
    # Black knights should be at A3, C3
    return (board[0][1] == 'B' and board[0][2] == 'B' and
            board[2][0] == 'w' and board[2][2] == 'w')

def find_solution():
    # Correct initial board state
    initial_board = [
        ['.', '.', '.', '.'],  # row 1
        ['.', '.', '.', '.'],  # row 2
        ['w', '.', 'w', '.']   # row 3
    ]
    initial_board[0][1] = 'B'  # B1
    initial_board[0][2] = 'B'  # C1
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {get_board_state(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        
        # Find all pieces of current color
        pieces = [(r, c) for r in range(3) for c in range(4) if board[r][c] == piece]
        
        for r, c in pieces:
            for new_r, new_c in get_knight_moves((r, c)):
                if board[new_r][new_c] == '.':
                    # Make move
                    new_board = [row[:] for row in board]
                    new_board[new_r][new_c] = piece
                    new_board[r][c] = '.'
                    
                    board_state = get_board_state(new_board)
                    if board_state not in visited:
                        visited.add(board_state)
                        move = f"{piece},{pos_to_coord((r,c))},{pos_to_coord((new_r,new_c))}"
                        queue.append((new_board, moves + [move], not is_white_turn))

    return None

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")