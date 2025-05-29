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
    return [(r,c) for r,c in moves if 0 <= r < 3 and 0 <= c < 3]

def pos_to_notation(pos):
    return f"{chr(65+pos[1])}{pos[0]+1}"

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # Check if white knights are at A1 and C1, black knights at C3 and A2
    return (board[0][0] == 'B' and board[0][2] == 'B' and 
            board[2][2] == 'w' and board[1][0] == 'w')

def find_solution():
    # Initial board
    initial_board = [
        ['B', '.', 'B'],
        ['w', '.', '.'],
        ['.', '.', 'w']
    ]
    
    # BFS queue: (board, moves_list, is_white_turn)
    queue = deque([(initial_board, [], True)])
    seen = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
            
        # Find current player's pieces
        piece = 'w' if is_white_turn else 'B'
        pieces = []
        for i in range(3):
            for j in range(3):
                if board[i][j] == piece:
                    pieces.append((i,j))
        
        # Try all possible moves for current pieces
        for start_pos in pieces:
            for end_pos in get_knight_moves(start_pos):
                if board[end_pos[0]][end_pos[1]] == '.':
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[end_pos[0]][end_pos[1]] = piece
                    new_board[start_pos[0]][start_pos[1]] = '.'
                    
                    board_str = board_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        move = f"{piece},{pos_to_notation(start_pos)},{pos_to_notation(end_pos)}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = find_solution()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")