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
    # White knights should be at A2, C1
    # Black knights should be at B2, D1
    white_target = {(0,1), (2,0)}  # A2, C1
    black_target = {(1,1), (3,0)}  # B2, D1
    
    current_white = set()
    current_black = set()
    
    for x in range(4):
        for y in range(3):
            if board[y][x] == 'w':
                current_white.add((x,y))
            elif board[y][x] == 'B':
                current_black.add((x,y))
    
    return current_white == white_target and current_black == black_target

def print_board(board):
    for row in board:
        print(row)
    print()

def find_solution():
    # Correct initial board state
    initial_board = [
        ['.', '.', '.', '.'],  # Row 3 (index 0)
        ['B', 'w', '.', '.'],  # Row 2 (index 1)
        ['.', '.', 'B', 'w']   # Row 1 (index 2)
    ]
    
    queue = deque([(initial_board, [], 'w')])
    visited = {board_to_string(initial_board)}
    
    while queue:
        board, moves, player = queue.popleft()
        
        if is_goal_state(board):
            return moves
        
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
                                from_sq = f"{chr(65+x)}{3-y}"  # Corrected coordinate system
                                to_sq = f"{chr(65+new_x)}{3-new_y}"
                                move = f"{player},{from_sq},{to_sq}"
                                next_player = 'B' if player == 'w' else 'w'
                                queue.append((new_board, moves + [move], next_player))
    
    return "No"

print(find_solution())