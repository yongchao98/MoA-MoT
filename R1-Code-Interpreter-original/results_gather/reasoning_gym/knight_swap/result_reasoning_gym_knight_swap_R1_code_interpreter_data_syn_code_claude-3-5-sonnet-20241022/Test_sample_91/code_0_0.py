from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1),
        (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2),
        (x-1, y+2), (x-1, y-2)
    ]
    return [(x, y) for x, y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # Check if white knights are at B3,D3 and black knights at A1,C1
    return (board[2][1] == 'w' and board[2][3] == 'w' and 
            board[0][0] == 'B' and board[0][2] == 'B')

def find_solution():
    # Initial board setup
    initial_board = [
        ['.', 'B', '.', 'B'],
        ['.', '.', '.', '.'],
        ['w', '.', 'w', '.']
    ]
    
    # BFS queue: (board, moves_list, is_white_turn)
    queue = deque([(initial_board, [], True)])
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
        
        # Find current player's knights
        piece = 'w' if is_white_turn else 'B'
        knights = []
        for y in range(3):
            for x in range(4):
                if board[y][x] == piece:
                    knights.append((x, y))
        
        # Try all possible moves for current player's knights
        for kx, ky in knights:
            for nx, ny in get_knight_moves((kx, ky)):
                if board[ny][nx] == '.':  # If target square is empty
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[ky][kx] = '.'
                    new_board[ny][nx] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        # Record move in chess notation
                        move = f"{piece},{chr(65+kx)}{ky+1},{chr(65+nx)}{ny+1}"
                        queue.append((new_board, moves + [move], not is_white_turn))

    return None

solution = find_solution()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")