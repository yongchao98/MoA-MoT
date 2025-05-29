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

def is_goal_state(board):
    # Check if white knights are at original black positions and vice versa
    return (board[2][1] == 'w' and board[0][0] == 'w' and 
            board[1][0] == 'B' and board[0][1] == 'B')

def find_knights(board):
    whites = []
    blacks = []
    for i in range(3):
        for j in range(4):
            if board[i][j] == 'w':
                whites.append((i,j))
            elif board[i][j] == 'B':
                blacks.append((i,j))
    return whites, blacks

def solve_knight_swap():
    initial_board = [
        ['.', 'B', '.', '.'],
        ['w', '.', '.', '.'],
        ['B', 'w', '.', '.']
    ]
    
    visited = set()
    queue = deque([(initial_board, [], True)])  # board, moves, white_turn
    
    while queue:
        board, moves, white_turn = queue.popleft()
        board_state = board_to_string(board)
        
        if board_state in visited:
            continue
            
        visited.add(board_state)
        
        if is_goal_state(board):
            return moves
            
        whites, blacks = find_knights(board)
        current_pieces = whites if white_turn else blacks
        piece_char = 'w' if white_turn else 'B'
        
        for piece_pos in current_pieces:
            r, c = piece_pos
            for new_r, new_c in get_knight_moves(piece_pos):
                if board[new_r][new_c] == '.':
                    new_board = copy.deepcopy(board)
                    new_board[r][c] = '.'
                    new_board[new_r][new_c] = piece_char
                    
                    move = f"{piece_char},{chr(c+65)}{3-r},{chr(new_c+65)}{3-new_r}"
                    queue.append((new_board, moves + [move], not white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")