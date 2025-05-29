from collections import deque
import copy

def get_knight_moves(pos):
    x, y = pos
    moves = [
        (x+2, y+1), (x+2, y-1), (x-2, y+1), (x-2, y-1),
        (x+1, y+2), (x+1, y-2), (x-1, y+2), (x-1, y-2)
    ]
    return [(x,y) for x,y in moves if 0 <= x < 4 and 0 <= y < 3]

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

def is_goal_state(board):
    # Check if white knights are at original black positions and vice versa
    black_positions = {(1,0), (2,0)}
    white_positions = {(1,1), (2,3)}
    
    current_black = set()
    current_white = set()
    for i in range(4):
        for j in range(3):
            if board[j][i] == 'B':
                current_black.add((i,j))
            elif board[j][i] == 'w':
                current_white.add((i,j))
    
    return current_white == black_positions and current_black == white_positions

def find_solution():
    # Initial board state
    initial_board = [
        ['.', '.', '.', '.'],
        ['B', 'w', '.', '.'],
        ['.', '.', 'B', 'w']
    ]
    
    queue = deque([(initial_board, [], 'w')])  # board, moves, next_player
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, player = queue.popleft()
        
        if is_goal_state(board):
            return moves
            
        # Find current player's pieces and possible moves
        for y in range(3):
            for x in range(4):
                if board[y][x] == player:
                    for new_x, new_y in get_knight_moves((x,y)):
                        if board[new_y][new_x] == '.':
                            new_board = copy.deepcopy(board)
                            new_board[y][x] = '.'
                            new_board[new_y][new_x] = player
                            board_str = board_to_string(new_board)
                            
                            if board_str not in visited:
                                visited.add(board_str)
                                move = f"{player},{chr(65+x)}{y+1},{chr(65+new_x)}{new_y+1}"
                                next_player = 'B' if player == 'w' else 'w'
                                queue.append((new_board, moves + [move], next_player))
    
    return None

solution = find_solution()
if solution:
    print(solution)
else:
    print("No")