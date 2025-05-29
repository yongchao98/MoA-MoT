from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and y >= 0 and y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_state(board):
    # White knights should be at C3, C2
    # Black knights should be at B1, D2
    return (board[0][2] == 'w' and 
            board[1][2] == 'w' and 
            board[2][1] == 'B' and 
            board[1][3] == 'B')

def validate_move(move, board, is_black_turn):
    piece = 'B' if is_black_turn else 'w'
    from_col = ord(move[2]) - ord('A')
    from_row = int(move[3]) - 1
    to_col = ord(move[5]) - ord('A')
    to_row = int(move[6]) - 1
    
    # Check if move is valid
    if not (0 <= from_col < 4 and 0 <= from_row < 3):
        return False
    if not (0 <= to_col < 4 and 0 <= to_row < 3):
        return False
    if board[from_row][from_col] != piece:
        return False
    if board[to_row][to_col] != '.':
        return False
    
    # Verify knight move
    dx = abs(to_col - from_col)
    dy = abs(to_row - from_row)
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def find_solution():
    initial_board = [
        ['.', '.', 'B', '.'],
        ['.', '.', 'B', 'w'],
        ['.', 'w', '.', '.']
    ]
    
    queue = deque([(initial_board, [], True)])
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_black_turn = queue.popleft()
        
        # Find current positions of pieces to move
        current_piece = 'B' if is_black_turn else 'w'
        for i in range(3):
            for j in range(4):
                if board[i][j] == current_piece:
                    for new_x, new_y in get_knight_moves((j, i)):
                        if board[new_y][new_x] == '.':
                            new_board = copy.deepcopy(board)
                            new_board[i][j] = '.'
                            new_board[new_y][new_x] = current_piece
                            
                            board_str = board_to_string(new_board)
                            if board_str not in visited:
                                visited.add(board_str)
                                move = f"{current_piece},{chr(65+j)}{i+1},{chr(65+new_x)}{new_y+1}"
                                new_moves = moves + [move]
                                
                                if is_target_state(new_board):
                                    # Validate entire sequence
                                    valid = True
                                    test_board = copy.deepcopy(initial_board)
                                    for m in new_moves:
                                        if not validate_move(m, test_board, m.startswith('B')):
                                            valid = False
                                            break
                                    if valid:
                                        return new_moves
                                
                                queue.append((new_board, new_moves, not is_black_turn))
    return "No"

solution = find_solution()
print(solution)