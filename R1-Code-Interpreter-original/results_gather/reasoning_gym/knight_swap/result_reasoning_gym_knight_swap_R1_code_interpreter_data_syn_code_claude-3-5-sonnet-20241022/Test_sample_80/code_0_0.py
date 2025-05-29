from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def pos_to_notation(pos):
    return chr(65 + pos[0]) + str(pos[1] + 1)

def state_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # Check if white knights are at original black positions and vice versa
    white_positions = set((i,j) for i in range(4) for j in range(3) 
                         if board[j][i] == 'w')
    black_positions = set((i,j) for i in range(4) for j in range(3) 
                         if board[j][i] == 'B')
    target_black = {(0,2), (1,0)}
    target_white = {(2,2), (3,1)}
    return white_positions == target_white and black_positions == target_black

def solve_knight_swap():
    # Initial board state
    initial_board = [
        ['w', '.', 'B', ' '],
        [' ', ' ', ' ', 'B'],
        [' ', 'w', ' ', '.']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    seen = {state_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        # Find all pieces of current color
        pieces = [(i,j) for i in range(4) for j in range(3) if board[j][i] == piece]
        
        for pos in pieces:
            for new_pos in get_knight_moves(pos):
                if board[new_pos[1]][new_pos[0]] == '.':
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[pos[1]][pos[0]] = '.'
                    new_board[new_pos[1]][new_pos[0]] = piece
                    
                    board_str = state_to_string(new_board)
                    if board_str not in seen:
                        seen.add(board_str)
                        move = f"{piece},{pos_to_notation(pos)},{pos_to_notation(new_pos)}"
                        queue.append((new_board, moves + [move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")