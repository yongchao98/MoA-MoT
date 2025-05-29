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
    return ''.join(''.join(row) for row in board)

def is_final_state(board):
    # Check if white knights are where black knights were and vice versa
    return (board[2][1] == 'w' and board[2][2] == 'w' and 
            board[1][1] == 'B' and board[0][3] == 'B')

def find_solution():
    # Initial board state
    initial_board = [
        ['.', 'B', 'B', ' '],
        ['.', 'w', '.', ' '],
        ['.', ' ', '.', 'w']
    ]
    
    # Queue for BFS: (board, moves_list, is_white_turn)
    queue = deque([(initial_board, [], True)])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, white_turn = queue.popleft()
        
        # Find current positions of pieces to move
        pieces = []
        for i in range(3):
            for j in range(4):
                if current_board[i][j] == ('w' if white_turn else 'B'):
                    pieces.append((i, j))
        
        # Try all possible moves for current pieces
        for piece in pieces:
            for new_pos in get_knight_moves(piece):
                if current_board[new_pos[0]][new_pos[1]] == '.':
                    # Make move
                    new_board = copy.deepcopy(current_board)
                    piece_type = new_board[piece[0]][piece[1]]
                    new_board[piece[0]][piece[1]] = '.'
                    new_board[new_pos[0]][new_pos[1]] = piece_type
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        move = f"{'w' if white_turn else 'B'},{chr(65+piece[1])}{piece[0]+1},{chr(65+new_pos[1])}{new_pos[0]+1}"
                        new_moves = moves + [move]
                        
                        if is_final_state(new_board):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not white_turn))
    
    return None

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")