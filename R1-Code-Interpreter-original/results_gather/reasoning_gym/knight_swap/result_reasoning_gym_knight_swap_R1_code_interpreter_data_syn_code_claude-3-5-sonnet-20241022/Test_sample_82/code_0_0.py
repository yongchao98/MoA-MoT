from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [(x+2,y+1), (x+2,y-1), (x-2,y+1), (x-2,y-1),
             (x+1,y+2), (x+1,y-2), (x-1,y+2), (x-1,y-2)]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 4]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_goal_state(board):
    # Check if white knights are at original black positions and vice versa
    return (board[0][0] == 'w' and board[0][1] == 'w' and 
            board[2][1] == 'B' and board[1][0] == 'B')

def find_solution():
    # Initial board state
    initial_board = [
        ['.', 'w', '.', '.'],
        ['w', '.', '.', '.'],
        ['B', 'B', '.', '.'],
        ['.', '.', '.', '.']
    ]
    
    # Queue stores (board, moves_list, is_white_turn)
    queue = deque([(initial_board, [], True)])
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, is_white_turn = queue.popleft()
        piece = 'w' if is_white_turn else 'B'
        
        # Find all pieces of current player
        pieces = []
        for i in range(4):
            for j in range(4):
                if board[i][j] == piece:
                    pieces.append((i,j))
        
        # Try all possible moves for each piece
        for pos in pieces:
            x, y = pos
            for new_x, new_y in get_knight_moves(pos):
                if board[new_x][new_y] == '.':  # If target square is empty
                    # Make move
                    new_board = copy.deepcopy(board)
                    new_board[x][y] = '.'
                    new_board[new_x][new_y] = piece
                    
                    board_str = board_to_string(new_board)
                    if board_str not in visited:
                        visited.add(board_str)
                        # Record move
                        move = f"{piece},{chr(y+65)}{x+1},{chr(new_y+65)}{new_x+1}"
                        new_moves = moves + [move]
                        
                        if is_goal_state(new_board):
                            return new_moves
                        
                        queue.append((new_board, new_moves, not is_white_turn))
    
    return None

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")