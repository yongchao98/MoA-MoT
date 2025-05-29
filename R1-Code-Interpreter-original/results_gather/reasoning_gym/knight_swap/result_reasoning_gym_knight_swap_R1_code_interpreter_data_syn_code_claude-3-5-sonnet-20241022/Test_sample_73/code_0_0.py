from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_goal_state(board):
    # Check if white knights are at B1, D1 and black knights at A2, C3
    return (board[0][2] == 'B' and board[1][0] == 'B' and 
            board[2][1] == 'w' and board[2][3] == 'w')

def find_solution():
    initial_board = [
        ['.', '.', 'w', '.'],
        ['w', '.', '.', '.'],
        ['.', 'B', '.', 'B']
    ]
    
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        piece = 'w' if is_white_turn else 'B'
        
        # Find all pieces of current color
        pieces = []
        for y in range(3):
            for x in range(4):
                if board[y][x] == piece:
                    pieces.append((x,y))
        
        # Try all possible moves for each piece
        for x,y in pieces:
            for new_x, new_y in get_knight_moves((x,y)):
                if board[new_y][new_x] == '.':
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[y][x] = '.'
                    new_board[new_y][new_x] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        
                        # Convert coordinates to chess notation
                        from_pos = f"{chr(65+x)}{y+1}"
                        to_pos = f"{chr(65+new_x)}{new_y+1}"
                        new_moves = moves + [f"{piece},{from_pos},{to_pos}"]
                        
                        if is_goal_state(new_board):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")