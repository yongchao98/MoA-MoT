from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_state(board):
    # Black knights should be at B1, D2
    # White knights should be at C2, C3
    return (board[2][1] == 'B' and  # B1
            board[1][3] == 'B' and  # D2
            board[1][2] == 'w' and  # C2
            board[0][2] == 'w')     # C3

def make_move(board, from_pos, to_pos, piece):
    new_board = copy.deepcopy(board)
    new_board[from_pos[1]][from_pos[0]] = '.'
    new_board[to_pos[1]][to_pos[0]] = piece
    return new_board

def find_solution():
    # Initial board state
    initial_board = [
        ['.', '.', 'B', '.'],  # Row 3 (index 0)
        ['.', '.', 'B', 'w'],  # Row 2 (index 1)
        ['.', 'w', '.', '.']   # Row 1 (index 2)
    ]
    
    # (board, moves_list, is_black_turn)
    queue = deque([(initial_board, [], True)])
    visited = {board_to_string(initial_board)}
    
    while queue:
        current_board, moves, is_black_turn = queue.popleft()
        current_piece = 'B' if is_black_turn else 'w'
        
        # Find all pieces of current color
        for y in range(3):
            for x in range(4):
                if current_board[y][x] == current_piece:
                    # Try all possible knight moves for this piece
                    for new_x, new_y in get_knight_moves((x, y)):
                        if current_board[new_y][new_x] == '.':
                            # Make the move
                            new_board = make_move(current_board, (x,y), (new_x,new_y), current_piece)
                            board_str = board_to_string(new_board)
                            
                            if board_str not in visited:
                                visited.add(board_str)
                                # Record move in algebraic notation
                                move = f"{current_piece},{chr(65+x)}{y+1},{chr(65+new_x)}{new_y+1}"
                                new_moves = moves + [move]
                                
                                if is_target_state(new_board):
                                    return new_moves
                                
                                queue.append((new_board, new_moves, not is_black_turn))
    
    return "No"

# Run the solution
solution = find_solution()
print(solution)