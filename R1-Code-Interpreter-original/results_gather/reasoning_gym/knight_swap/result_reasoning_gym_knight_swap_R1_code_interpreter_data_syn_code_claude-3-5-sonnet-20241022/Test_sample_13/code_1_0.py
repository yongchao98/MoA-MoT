from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def pos_to_notation(pos):
    return chr(65 + pos[0]) + str(3 - pos[1])  # Fixed coordinate system

def notation_to_pos(notation):
    return (ord(notation[0]) - 65, 3 - int(notation[1]))

def is_target_reached(board):
    # Check if white knights are at B1,C1 and black knights are at D2,D3
    return (board[2][1] == 'w' and board[2][2] == 'w' and 
            board[0][3] == 'B' and board[1][3] == 'B')

def get_board_state(board):
    return ''.join(''.join(row) for row in board)

def solve_knight_swap():
    # Initial board setup (0-based indexing, reversed rows)
    initial_board = [
        ['.','.','.' ,'w'],  # row 3 (index 0)
        ['.','.','.', 'w'],  # row 2 (index 1)
        ['.','B','B', '.']   # row 1 (index 2)
    ]
    
    queue = deque([(initial_board, [], True)])  # (board, moves, is_white_turn)
    visited = {get_board_state(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        
        # Find all pieces of current color
        pieces = []
        for y in range(3):
            for x in range(4):
                if board[y][x] == piece:
                    pieces.append((x,y))
        
        # Try all possible moves for each piece
        for piece_pos in pieces:
            x, y = piece_pos
            for new_x, new_y in get_knight_moves(piece_pos):
                if board[new_y][new_x] == '.':  # If target square is empty
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[y][x] = '.'
                    new_board[new_y][new_x] = piece
                    
                    board_state = get_board_state(new_board)
                    if board_state not in visited:
                        visited.add(board_state)
                        move = f"{piece},{pos_to_notation((x,y))},{pos_to_notation((new_x,new_y))}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")