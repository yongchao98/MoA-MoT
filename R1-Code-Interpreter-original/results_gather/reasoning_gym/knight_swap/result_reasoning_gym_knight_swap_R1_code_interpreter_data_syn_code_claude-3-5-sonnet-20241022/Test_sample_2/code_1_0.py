from collections import deque

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

def is_target_reached(board):
    # White knights should be at A2, B1 (original black positions)
    # Black knights should be at C3, D1 (original white positions)
    return (board[1][0] == 'w' and board[0][1] == 'w' and 
            board[2][2] == 'B' and board[2][3] == 'B')

def get_valid_moves(board, is_white_turn):
    piece = 'w' if is_white_turn else 'B'
    moves = []
    for r in range(3):
        for c in range(4):
            if board[r][c] == piece:
                for new_r, new_c in get_knight_moves((r,c)):
                    if board[new_r][new_c] == '.':
                        moves.append(((r,c), (new_r,new_c)))
    return moves

def pos_to_notation(pos):
    return chr(pos[1] + ord('A')) + str(pos[0] + 1)

def print_board(board):
    for row in board:
        print(row)
    print()

def solve_knight_swap():
    # Initial board state
    initial_board = [
        ['.', '.', 'w', '.'],
        ['B', '.', '.', '.'],
        ['.', 'B', '.', 'w']
    ]
    
    queue = deque([(initial_board, [], True)])
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if len(moves) > 15:  # Limit search depth
            continue
            
        if is_target_reached(board):
            return moves
            
        valid_moves = get_valid_moves(board, is_white_turn)
        for (r,c), (new_r,new_c) in valid_moves:
            new_board = [row[:] for row in board]
            new_board[new_r][new_c] = new_board[r][c]
            new_board[r][c] = '.'
            
            board_str = board_to_string(new_board)
            if board_str not in visited:
                visited.add(board_str)
                piece = 'w' if is_white_turn else 'B'
                move_notation = f"{piece},{pos_to_notation((r,c))},{pos_to_notation((new_r,new_c))}"
                queue.append((new_board, moves + [move_notation], not is_white_turn))
    
    return "No"

solution = solve_knight_swap()
print(solution)