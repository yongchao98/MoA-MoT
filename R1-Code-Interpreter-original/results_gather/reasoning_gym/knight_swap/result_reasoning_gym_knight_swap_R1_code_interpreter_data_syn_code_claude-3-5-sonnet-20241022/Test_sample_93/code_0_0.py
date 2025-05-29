from collections import deque
import copy

def get_knight_moves(pos):
    row, col = pos
    moves = [
        (row+2, col+1), (row+2, col-1),
        (row-2, col+1), (row-2, col-1),
        (row+1, col+2), (row+1, col-2),
        (row-1, col+2), (row-1, col-2)
    ]
    return [(r,c) for r,c in moves if 0 <= r < 3 and 0 <= c < 4]

def board_to_string(board):
    return ','.join(''.join(row) for row in board)

def is_target_reached(board):
    # Check if white knights are at original black positions and vice versa
    return (board[2][1] == 'w' and board[0][0] == 'w' and 
            board[1][0] == 'B' and board[0][1] == 'B')

def find_solution():
    # Initial board
    initial_board = [
        ['.', 'B', '.', '.'],
        ['w', '.', '.', '.'],
        ['B', 'w', '.', '.']
    ]
    
    # BFS queue: (board, moves_list, is_white_turn)
    queue = deque([(initial_board, [], True)])
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        
        # Find all pieces of current color
        pieces = []
        for r in range(3):
            for c in range(4):
                if board[r][c] == piece:
                    pieces.append((r,c))
        
        # Try all possible moves for each piece
        for r,c in pieces:
            for new_r, new_c in get_knight_moves((r,c)):
                if board[new_r][new_c] == '.':
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[r][c] = '.'
                    new_board[new_r][new_c] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        # Convert coordinates to chess notation
                        from_pos = f"{chr(65+c)}{3-r}"
                        to_pos = f"{chr(65+new_c)}{3-new_r}"
                        new_moves = moves + [f"{piece},{from_pos},{to_pos}"]
                        queue.append((new_board, new_moves, not is_white_turn))

    return None

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")