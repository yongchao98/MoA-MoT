from collections import deque
import copy

def is_valid(x, y, rows, cols, board):
    return 0 <= x < rows and 0 <= y < cols and board[x][y] != '+'

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def is_box(char):
    return char in ['@', '$']

def is_goal(char):
    return char in ['X', '%', '$']

def is_complete(board):
    goals = 0
    boxes_on_goals = 0
    for row in board:
        for cell in row:
            if is_goal(cell):
                goals += 1
            if cell == '$':
                boxes_on_goals += 1
    return goals == boxes_on_goals

def board_to_string(board):
    return '\n'.join([''.join(row) for row in board])

def solve_sokoban(initial_board):
    moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    queue = deque([(initial_board, "")])
    visited = set()
    
    while queue:
        board, path = queue.popleft()
        board_str = board_to_string(board)
        
        if board_str in visited:
            continue
            
        visited.add(board_str)
        
        if is_complete(board):
            return path
            
        px, py = get_player_pos(board)
        
        for move, (dx, dy) in moves.items():
            new_px, new_py = px + dx, py + dy
            
            if not is_valid(new_px, new_py, len(board), len(board[0]), board):
                continue
                
            new_board = [list(row) for row in board]
            
            if is_box(board[new_px][new_py]):
                box_px, box_py = new_px + dx, new_py + dy
                
                if not is_valid(box_px, box_py, len(board), len(board[0]), board) or is_box(board[box_px][box_py]):
                    continue
                    
                # Move box
                new_board[box_px][box_py] = '$' if is_goal(board[box_px][box_py]) else '@'
                
            # Move player
            new_board[new_px][new_py] = '%' if is_goal(board[new_px][new_py]) else '*'
            new_board[px][py] = 'X' if is_goal(board[px][py]) else '-'
            
            queue.append((new_board, path + move))
    
    return None

# Initial board
initial_board = [
    list("+++++++++++"),
    list("+-------+++"),
    list("+--@----+++"),
    list("+---XX@--+"),
    list("+-$$$-@-@+"),
    list("+------X%+"),
    list("+++++++++++")
]

solution = solve_sokoban(initial_board)
print(solution)