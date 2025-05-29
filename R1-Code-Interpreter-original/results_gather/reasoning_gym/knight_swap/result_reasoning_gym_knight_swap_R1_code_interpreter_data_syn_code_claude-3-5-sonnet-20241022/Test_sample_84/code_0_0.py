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
    # Check if white knights are at black's initial positions and vice versa
    return (board[0][2] == 'w' and board[1][3] == 'w' and 
            board[1][2] == 'B' and board[1][1] == 'B')

def find_solution():
    # Initial board state
    initial_board = [
        ['.', '.', 'B', '.'],
        ['.', '.', 'w', 'B'],
        ['.', 'w', '.', '.']
    ]
    
    # Queue stores (board, moves_list, is_white_turn)
    queue = deque([(initial_board, [], True)])
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, white_turn = queue.popleft()
        
        # Find current player's pieces
        piece = 'w' if white_turn else 'B'
        pieces = []
        for i in range(3):
            for j in range(4):
                if board[i][j] == piece:
                    pieces.append((j, i))  # (x,y) format
        
        # Try all possible moves for each piece
        for piece_pos in pieces:
            x, y = piece_pos
            for new_x, new_y in get_knight_moves(piece_pos):
                if board[new_y][new_x] == '.':  # If target square is empty
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[y][x] = '.'
                    new_board[new_y][new_x] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        # Record move
                        from_pos = f"{chr(65+x)}{3-y}"
                        to_pos = f"{chr(65+new_x)}{3-new_y}"
                        new_moves = moves + [f"{piece},{from_pos},{to_pos}"]
                        
                        if is_goal_state(new_board):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not white_turn))
    
    return None

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")