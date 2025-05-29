from collections import deque
import copy

def get_knight_moves(x, y):
    moves = [
        (x+2, y+1), (x+2, y-1),
        (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2),
        (x-1, y+2), (x-1, y-2)
    ]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_target_state(board):
    # Check if white knights are at A3 and B1, black knights at C2 and D1
    return (board[0][0] == 'w' and board[2][1] == 'w' and 
            board[1][2] == 'B' and board[2][3] == 'B')

def solve_knight_swap():
    # Initial board state (0,0 is A1)
    initial_board = [
        ['.','.','w','.'],  # Row 1
        ['.','.','.','w'],  # Row 2
        ['B','B','.','.']]  # Row 3
    
    queue = deque([(initial_board, [], True)])  # board, moves, white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, white_turn = queue.popleft()
        
        if is_target_state(board):
            return moves
        
        piece = 'w' if white_turn else 'B'
        
        # Find all pieces of current color
        for i in range(3):
            for j in range(4):
                if board[i][j] == piece:
                    # Try all possible moves for this piece
                    for new_x, new_y in get_knight_moves(i, j):
                        if board[new_x][new_y] == '.':
                            # Make move
                            new_board = copy.deepcopy(board)
                            new_board[i][j] = '.'
                            new_board[new_x][new_y] = piece
                            
                            board_str = board_to_string(new_board)
                            if board_str not in visited:
                                visited.add(board_str)
                                # Convert coordinates to chess notation
                                from_pos = f"{chr(65+j)}{3-i}"
                                to_pos = f"{chr(65+new_y)}{3-new_x}"
                                new_moves = moves + [f"{piece},{from_pos},{to_pos}"]
                                queue.append((new_board, new_moves, not white_turn))
    
    return None

solution = solve_knight_swap()
if solution:
    print(solution)
else:
    print("No")