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
    return ','.join(''.join(row) for row in board)

def is_target_reached(board):
    # Check if white knights are at B3 and D2, black knights at A1 and C1
    return (board[0][1] == 'w' and board[1][3] == 'w' and 
            board[2][0] == 'B' and board[2][2] == 'B')

def find_solution():
    # Initial board state (0,0 is A3)
    initial_board = [
        ['.', 'B', '.', '.'],
        ['.', '.', '.', 'B'],
        ['w', '.', 'w', '.']
    ]
    
    # Using BFS to find solution
    queue = deque([(initial_board, [], True)])  # board, moves, is_white_turn
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        
        if is_target_reached(board):
            return moves
            
        piece = 'w' if is_white_turn else 'B'
        
        # Find all pieces of current color
        pieces = []
        for i in range(3):
            for j in range(4):
                if board[i][j] == piece:
                    pieces.append((j, i))  # x, y coordinates
        
        # Try all possible moves for each piece
        for piece_pos in pieces:
            for new_pos in get_knight_moves(piece_pos):
                if board[new_pos[1]][new_pos[0]] == '.':
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[piece_pos[1]][piece_pos[0]] = '.'
                    new_board[new_pos[1]][new_pos[0]] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        # Convert coordinates to chess notation
                        from_sq = f"{chr(65+piece_pos[0])}{3-piece_pos[1]}"
                        to_sq = f"{chr(65+new_pos[0])}{3-new_pos[1]}"
                        new_moves = moves + [f"{piece},{from_sq},{to_sq}"]
                        queue.append((new_board, new_moves, not is_white_turn))

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")