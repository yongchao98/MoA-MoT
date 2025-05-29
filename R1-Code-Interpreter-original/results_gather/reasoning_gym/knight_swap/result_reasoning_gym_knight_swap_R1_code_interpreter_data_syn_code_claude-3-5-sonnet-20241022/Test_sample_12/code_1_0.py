from collections import deque
import copy

def get_knight_moves(pos):
    moves = [
        (2,1), (2,-1), (-2,1), (-2,-1),
        (1,2), (1,-2), (-1,2), (-1,-2)
    ]
    x, y = pos
    valid = []
    for dx, dy in moves:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 3 and 0 <= new_y < 3:
            valid.append((new_x, new_y))
    return valid

def pos_to_coord(pos):
    return f"{chr(65+pos[1])}{3-pos[0]}"  # Fixed coordinate conversion

def coord_to_pos(coord):
    return (3-int(coord[1]), ord(coord[0])-65)

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_state(board):
    # Check if white knights are at C1 and C3, black knights at A2 and B1
    return (board[2][2] == 'w' and board[0][2] == 'w' and 
            board[1][0] == 'B' and board[2][1] == 'B')

def find_solution():
    # Initial board state (0,0 is top-left)
    initial = [
        ['.', ' ', 'B'],  # row 3
        ['w', ' ', '.'],  # row 2
        [' ', 'w', 'B']   # row 1
    ]
    
    queue = deque([(initial, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_state(board):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        
        # Find all pieces of current color
        pieces = []
        for i in range(3):
            for j in range(3):
                if board[i][j] == piece:
                    pieces.append((i,j))
        
        # Try all possible moves for each piece
        for piece_pos in pieces:
            for new_pos in get_knight_moves(piece_pos):
                i, j = new_pos
                if board[i][j] in ['.', ' ']:
                    new_board = copy.deepcopy(board)
                    # Make the move
                    new_board[i][j] = piece
                    new_board[piece_pos[0]][piece_pos[1]] = '.' if (piece_pos[0] + piece_pos[1]) % 2 == 0 else ' '
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        move = f"{piece},{pos_to_coord(piece_pos)},{pos_to_coord(new_pos)}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = find_solution()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")