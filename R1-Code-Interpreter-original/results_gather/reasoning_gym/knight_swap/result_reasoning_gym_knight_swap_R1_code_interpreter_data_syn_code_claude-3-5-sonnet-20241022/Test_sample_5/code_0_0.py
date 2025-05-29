from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def pos_to_notation(pos):
    return chr(65 + pos[0]) + str(pos[1] + 1)

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_reached(board):
    # Check if black knights are at white starting positions and vice versa
    return (board[0][0] == 'B' and board[0][2] == 'B' and 
            board[1][1] == 'w' and board[1][3] == 'w')

def solve_knight_swap():
    # Initial board state
    initial_board = [
        ['.', '.', '.', '.'],
        ['.', 'B', '.', 'B'],
        ['w', '.', 'w', '.']
    ]
    
    seen = set()
    queue = deque([(initial_board, [], True)])  # board, moves, is_black_turn
    
    while queue:
        board, moves, is_black_turn = queue.popleft()
        board_str = board_to_string(board)
        
        if board_str in seen:
            continue
        seen.add(board_str)
        
        if is_target_reached(board):
            return moves
            
        piece = 'B' if is_black_turn else 'w'
        
        # Find all pieces of current color
        pieces = [(r,c) for r in range(3) for c in range(4) 
                 if board[r][c] == piece]
        
        # Try all possible moves for each piece
        for pos in pieces:
            for new_pos in get_knight_moves(pos):
                if board[new_pos[0]][new_pos[1]] == '.':
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[pos[0]][pos[1]] = '.'
                    new_board[new_pos[0]][new_pos[1]] = piece
                    
                    # Record move
                    move = f"{piece},{pos_to_notation(pos)},{pos_to_notation(new_pos)}"
                    new_moves = moves + [move]
                    
                    queue.append((new_board, new_moves, not is_black_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")