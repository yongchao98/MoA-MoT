from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(x, y) for x, y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_goal_state(board):
    # Check if white knights are at B3 and C1, black knights at B1 and D2
    return (board[0][1] == 'w' and board[2][2] == 'w' and 
            board[2][1] == 'B' and board[1][3] == 'B')

def get_piece_positions(board, piece):
    return [(i, j) for i in range(4) for j in range(3) 
            if board[j][i] == piece]

def solve_knight_swap():
    # Initial board state
    initial_board = [
        ['.', 'B', '.', ''],
        ['.', '.', '.', 'w'],
        ['.', 'w', 'B', '']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_goal_state(board):
            return moves
            
        current_piece = 'w' if is_white_turn else 'B'
        pieces = get_piece_positions(board, current_piece)
        
        for piece_pos in pieces:
            valid_moves = get_knight_moves(piece_pos)
            for new_pos in valid_moves:
                if board[new_pos[1]][new_pos[0]] == '.':
                    new_board = [row[:] for row in board]
                    # Make the move
                    new_board[piece_pos[1]][piece_pos[0]] = '.'
                    new_board[new_pos[1]][new_pos[0]] = current_piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        move = f"{current_piece},{chr(65+piece_pos[0])}{piece_pos[1]+1},{chr(65+new_pos[0])}{new_pos[1]+1}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")