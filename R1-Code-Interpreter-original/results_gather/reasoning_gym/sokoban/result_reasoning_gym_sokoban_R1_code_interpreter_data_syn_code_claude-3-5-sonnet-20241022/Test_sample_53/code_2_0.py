from collections import deque

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                return i, j
    return None

def is_goal_position(board, x, y):
    return board[x][y] in ['X', '%', '$']

def make_move(board, player_x, player_y, dx, dy):
    new_board = [row[:] for row in board]
    new_x, new_y = player_x + dx, player_y + dy
    
    # Check if the move is to a wall
    if new_board[new_x][new_y] == '+':
        return None

    # Moving to empty space or goal
    if new_board[new_x][new_y] in ['-', 'X']:
        was_on_goal = is_goal_position(board, player_x, player_y)
        new_board[player_x][player_y] = 'X' if was_on_goal else '-'
        new_board[new_x][new_y] = '%' if new_board[new_x][new_y] == 'X' else '*'
        return new_board

    # Moving to a box
    if new_board[new_x][new_y] in ['@', '$']:
        box_x, box_y = new_x + dx, new_y + dy
        
        # Check if box can be pushed
        if not is_valid(box_x, box_y, len(board), len(board[0])) or \
           new_board[box_x][box_y] in ['+', '@', '$']:
            return None
            
        was_on_goal = is_goal_position(board, player_x, player_y)
        box_to_goal = is_goal_position(board, box_x, box_y)
        
        new_board[player_x][player_y] = 'X' if was_on_goal else '-'
        new_board[new_x][new_y] = '%' if is_goal_position(board, new_x, new_y) else '*'
        new_board[box_x][box_y] = '$' if box_to_goal else '@'
        
        return new_board
        
    return None

def is_solved(board):
    goals = 0
    boxes_on_goals = 0
    for row in board:
        for cell in row:
            if cell in ['X', '%']:
                goals += 1
            if cell == '$':
                boxes_on_goals += 1
    return goals == 0 and boxes_on_goals > 0

def solve_sokoban(board):
    moves = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    queue = deque([(board, "")])
    seen = {str(board)}
    
    while queue:
        current_board, path = queue.popleft()
        
        if is_solved(current_board):
            return path
            
        player_x, player_y = get_player_pos(current_board)
        
        for move_char, dx, dy in moves:
            new_board = make_move(current_board, player_x, player_y, dx, dy)
            
            if new_board is not None:
                board_str = str(new_board)
                if board_str not in seen:
                    seen.add(board_str)
                    queue.append((new_board, path + move_char))
    
    return None

# Initialize the board
initial_board = [
    list("+++++++" ),
    list("+++++X+" ),
    list("++++-@+" ),
    list("+++-@*+" ),
    list("++$-X$+" ),
    list("+-+-+++" ),
    list("+@+++++" ),
    list("+X+++++" ),
    list("+++$+$+" ),
    list("+++++++" )
]

solution = solve_sokoban(initial_board)
print(solution if solution else "No solution found")