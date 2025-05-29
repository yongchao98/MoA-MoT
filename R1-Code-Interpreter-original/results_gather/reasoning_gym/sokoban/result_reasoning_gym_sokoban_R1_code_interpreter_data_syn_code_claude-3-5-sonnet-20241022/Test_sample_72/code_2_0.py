def get_player_pos(board):
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] in ['*', '%']:
                return (i, j)
    return None

def is_goal_state(board):
    boxes = 0
    goals = 0
    for row in board:
        for cell in row:
            if cell in ['@', '$']:
                boxes += 1
            if cell in ['X', '$']:
                goals += 1
    return boxes == 0 and goals > 0

def solve_sokoban(board, max_moves=25):
    def try_move(board, pos, direction):
        px, py = pos
        dx, dy = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}[direction]
        nx, ny = px + dx, py + dy
        
        # Copy board
        new_board = [row[:] for row in board]
        
        # Check if move is valid
        if new_board[nx][ny] == '+':
            return None
            
        # Moving to empty space or goal
        if new_board[nx][ny] in ['-', 'X']:
            new_board[px][py] = 'X' if board[px][py] == '%' else '-'
            new_board[nx][ny] = '%' if new_board[nx][ny] == 'X' else '*'
            return new_board
            
        # Pushing box
        if new_board[nx][ny] in ['@', '$']:
            nnx, nny = nx + dx, ny + dy
            if (0 <= nnx < len(board) and 0 <= nny < len(board[0]) and 
                new_board[nnx][nny] in ['-', 'X']):
                new_board[px][py] = 'X' if board[px][py] == '%' else '-'
                new_board[nx][ny] = '%' if new_board[nx][ny] == '$' else '*'
                new_board[nnx][nny] = '$' if new_board[nnx][nny] == 'X' else '@'
                return new_board
        return None

    def dfs(board, moves, path):
        if len(path) > max_moves:
            return None
            
        if is_goal_state(board):
            return path
            
        player_pos = get_player_pos(board)
        if not player_pos:
            return None
            
        for direction in ['U', 'D', 'L', 'R']:
            new_board = try_move(board, player_pos, direction)
            if new_board:
                board_str = str(new_board)
                if board_str not in moves:
                    moves.add(board_str)
                    result = dfs(new_board, moves, path + direction)
                    if result:
                        return result
                    moves.remove(board_str)
        return None

    initial_moves = {str(board)}
    return dfs(board, initial_moves, "")

# Initial board
initial_board = [
    ['+', '+', '+', '+', '+', '+', '+', '+', '+'],
    ['+', '-', '-', '-', '-', '-', '-', 'X', '+'],
    ['+', '@', '-', '-', '-', 'X', '$', '-', '+'],
    ['+', '%', '@', '-', '-', '-', '@', '-', '+'],
    ['+', '@', '-', '-', '-', '@', 'X', '-', '+'],
    ['+', '-', '@', '-', '-', '-', 'X', 'X', '+'],
    ['+', 'X', '+', '-', '-', '-', '-', '@', '+'],
    ['+', '-', '-', '-', '-', '-', '-', '-', '+'],
    ['+', '+', '+', '+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_board)
print(f"<<<{solution}>>>")