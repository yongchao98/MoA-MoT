from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0<=x<4 and 0<=y<3]

def pos_to_notation(pos):
    return chr(pos[0] + ord('A')) + str(pos[1] + 1)

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_state(board):
    # Check if white knights are at original black positions and vice versa
    black_positions = {(0,2), (0,1)}  # A3, A2
    white_positions = {(1,0), (3,1)}  # B1, D2
    
    for i in range(4):
        for j in range(3):
            if board[j][i] == 'w' and (i,j) not in black_positions:
                return False
            if board[j][i] == 'B' and (i,j) not in white_positions:
                return False
    return True

def solve_puzzle():
    # Initial board state
    initial_board = [
        ['B', '.', '.', '.'],
        ['B', '.', '.', 'w'],
        ['.', 'w', '.', '.']
    ]
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        board_str = board_to_string(board)
        
        if board_str in visited:
            continue
            
        visited.add(board_str)
        
        if is_target_state(board):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        
        # Find all pieces of current player
        for y in range(3):
            for x in range(4):
                if board[y][x] == piece:
                    # Try all possible moves for this piece
                    for new_x, new_y in get_knight_moves((x,y)):
                        if board[new_y][new_x] == '.':
                            # Make move
                            new_board = copy.deepcopy(board)
                            new_board[y][x] = '.'
                            new_board[new_y][new_x] = piece
                            
                            move = f"{piece},{pos_to_notation((x,y))},{pos_to_notation((new_x,new_y))}"
                            queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_puzzle()
if solution:
    print(solution)
else:
    print("No")