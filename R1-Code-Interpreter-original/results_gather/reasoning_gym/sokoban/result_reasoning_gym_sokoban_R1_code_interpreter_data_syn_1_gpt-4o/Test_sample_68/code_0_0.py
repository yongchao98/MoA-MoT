from collections import deque

# Define the initial board state
initial_board = [
    "+ + + + + + + + +",
    "+ - - - - X + + +",
    "+ - - - - @ X + +",
    "+ + - - @ X $ + +",
    "+ + - + - - + + +",
    "+ + - $ + @ $ X +",
    "+ + + + - - + @ +",
    "+ + + X - - @ * +",
    "+ + + + + + + + +"
]

# Convert the board to a more manageable format
board = [list(row.replace(" ", "")) for row in initial_board]

# Define the directions for movement
directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}

# Find the initial position of the player
def find_player(board):
    for r, row in enumerate(board):
        for c, cell in enumerate(row):
            if cell == '*':
                return r, c
    return None

# Check if a position is within the board and not a wall
def is_valid_position(r, c, board):
    return 0 <= r < len(board) and 0 <= c < len(board[0]) and board[r][c] != '+'

# Check if the puzzle is solved
def is_solved(board):
    for row in board:
        for cell in row:
            if cell == '@':
                return False
    return True

# Move the player and update the board
def move_player(r, c, dr, dc, board):
    new_r, new_c = r + dr, c + dc
    if not is_valid_position(new_r, new_c, board):
        return None, None
    if board[new_r][new_c] in '-X':
        board[r][c], board[new_r][new_c] = '-', '*'
        return new_r, new_c
    if board[new_r][new_c] in '@$':
        box_r, box_c = new_r + dr, new_c + dc
        if not is_valid_position(box_r, box_c, board) or board[box_r][box_c] not in '-X':
            return None, None
        board[r][c], board[new_r][new_c], board[box_r][box_c] = '-', '*', '@' if board[box_r][box_c] == '-' else '$'
        return new_r, new_c
    return None, None

# Perform BFS to find the solution
def solve_sokoban(board):
    start_r, start_c = find_player(board)
    queue = deque([(start_r, start_c, board, "")])
    visited = set()
    visited.add(tuple(map(tuple, board)))
    
    while queue:
        r, c, current_board, path = queue.popleft()
        
        if is_solved(current_board):
            return path
        
        for move, (dr, dc) in directions.items():
            new_board = [row[:] for row in current_board]
            new_r, new_c = move_player(r, c, dr, dc, new_board)
            if new_r is not None and new_c is not None:
                board_tuple = tuple(map(tuple, new_board))
                if board_tuple not in visited:
                    visited.add(board_tuple)
                    queue.append((new_r, new_c, new_board, path + move))
    
    return "No solution"

# Find the solution
solution = solve_sokoban(board)
print(solution)