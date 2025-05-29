from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def pos_to_algebraic(pos):
    return chr(pos[0] + ord('A')) + str(pos[1] + 1)

def state_to_string(board):
    return ''.join(''.join(row) for row in board)

def find_solution():
    # Initial board state
    initial_board = [
        ['.', '.', 'B', '.'],
        ['w', 'B', '.', '.'],
        ['.', '.', '.', 'w']
    ]
    
    # Target positions for verification
    white_targets = {(1,1), (2,0)}  # B2, C3 positions
    black_targets = {(0,1), (3,0)}  # A2, D1 positions
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    seen = {state_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        # Check if we reached target state
        whites = {(x,y) for y in range(3) for x in range(4) if board[y][x] == 'w'}
        blacks = {(x,y) for y in range(3) for x in range(4) if board[y][x] == 'B'}
        
        if whites == white_targets and blacks == black_targets:
            return moves
            
        # Generate next moves
        pieces = whites if is_white_turn else blacks
        piece_char = 'w' if is_white_turn else 'B'
        
        for piece_pos in pieces:
            for new_pos in get_knight_moves(piece_pos):
                if board[new_pos[1]][new_pos[0]] == '.':
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[piece_pos[1]][piece_pos[0]] = '.'
                    new_board[new_pos[1]][new_pos[0]] = piece_char
                    
                    board_str = state_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        move = f"{piece_char},{pos_to_algebraic(piece_pos)},{pos_to_algebraic(new_pos)}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")