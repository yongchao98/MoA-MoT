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

def is_target_state(board):
    # Check if white knights are at original black positions and vice versa
    return (board[0][1] == 'w' and board[0][2] == 'w' and 
            board[1][1] == 'B' and board[0][3] == 'B')

def get_valid_moves(board, is_white_turn):
    piece = 'w' if is_white_turn else 'B'
    moves = []
    
    # Find all pieces of current color
    for i in range(3):
        for j in range(4):
            if board[i][j] == piece:
                for new_r, new_c in get_knight_moves((i,j)):
                    if board[new_r][new_c] == '.':
                        moves.append(((i,j), (new_r,new_c)))
    return moves

def pos_to_notation(pos):
    return chr(pos[1] + ord('A')) + str(3 - pos[0])

def solve_knight_swap():
    initial_board = [
        ['.', 'B', 'B', 'w'],
        ['.', 'w', '.', '.'],
        ['.', '.', '.', '.']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_state(board):
            return moves
            
        valid_moves = get_valid_moves(board, is_white_turn)
        
        for (from_r, from_c), (to_r, to_c) in valid_moves:
            new_board = [row[:] for row in board]
            new_board[to_r][to_c] = new_board[from_r][from_c]
            new_board[from_r][from_c] = '.'
            
            board_str = board_to_string(new_board)
            if board_str not in visited:
                visited.add(board_str)
                piece = 'w' if is_white_turn else 'B'
                move_notation = f"{piece},{pos_to_notation((from_r,from_c))},{pos_to_notation((to_r,to_c))}"
                queue.append((new_board, moves + [move_notation], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")