from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_valid_move(board, from_pos, to_pos, piece):
    # Check if the move is within bounds
    if not (0 <= from_pos[0] < 4 and 0 <= from_pos[1] < 3 and
            0 <= to_pos[0] < 4 and 0 <= to_pos[1] < 3):
        return False
    
    # Check if the piece is at the starting position
    if board[from_pos[1]][from_pos[0]] != piece:
        return False
    
    # Check if the destination is empty
    if board[to_pos[1]][to_pos[0]] != '.':
        return False
    
    # Check if it's a valid knight move
    dx = abs(to_pos[0] - from_pos[0])
    dy = abs(to_pos[1] - from_pos[1])
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def is_target_state(board):
    # Black knights should be at B1 and D2
    # White knights should be at C2 and C3
    return (board[2][1] == 'B' and  # B1
            board[1][3] == 'B' and  # D2
            board[1][2] == 'w' and  # C2
            board[0][2] == 'w')     # C3

def find_solution():
    # Initial board state (0,0 is top-left)
    initial_board = [
        ['.', '.', 'B', '.'],  # Row 3 (index 0)
        ['.', '.', 'B', 'w'],  # Row 2 (index 1)
        ['.', 'w', '.', '.']   # Row 1 (index 2)
    ]
    
    queue = deque([(initial_board, [], True)])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, is_black_turn = queue.popleft()
        piece = 'B' if is_black_turn else 'w'
        
        # Find all pieces of current color
        for y in range(3):
            for x in range(4):
                if current_board[y][x] == piece:
                    # Try all possible knight moves
                    for new_x, new_y in get_knight_moves((x, y)):
                        if is_valid_move(current_board, (x,y), (new_x,new_y), piece):
                            # Make move
                            new_board = copy.deepcopy(current_board)
                            new_board[y][x] = '.'
                            new_board[new_y][new_x] = piece
                            
                            board_str = board_to_string(new_board)
                            if board_str not in visited:
                                visited.add(board_str)
                                move = f"{piece},{chr(65+x)}{3-y},{chr(65+new_x)}{3-new_y}"
                                new_moves = moves + [move]
                                
                                if is_target_state(new_board):
                                    return new_moves
                                
                                queue.append((new_board, new_moves, not is_black_turn))
    return "No"

solution = find_solution()
print(solution)