def is_valid_pos(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_pos(x, y, direction):
    if direction == 'U': return (x-1, y)
    if direction == 'D': return (x+1, y)
    if direction == 'L': return (x, y-1)
    if direction == 'R': return (x, y+1)

def try_move(board, player_pos, direction, depth=0, max_depth=50, visited=None):
    if depth >= max_depth:
        return None
    
    if visited is None:
        visited = set()
    
    rows, cols = len(board), len(board[0])
    x, y = player_pos
    new_x, new_y = get_next_pos(x, y, direction)
    
    if not is_valid_pos(new_x, new_y, rows, cols) or board[new_x][new_y] == '+':
        return None
    
    new_board = [list(row) for row in board]
    
    # If moving to empty space or goal
    if board[new_x][new_y] in ['-', 'X']:
        new_board[x][y] = '-' if board[x][y] == '*' else 'X'
        new_board[new_x][new_y] = '*'
        state = (''.join(''.join(row) for row in new_board))
        if state in visited:
            return None
        visited.add(state)
        return (new_board, (new_x, new_y))
    
    # If pushing a box
    if board[new_x][new_y] in ['@', '$']:
        box_x, box_y = get_next_pos(new_x, new_y, direction)
        if not is_valid_pos(box_x, box_y, rows, cols) or board[box_x][box_y] in ['+', '@', '$']:
            return None
            
        new_board[x][y] = '-' if board[x][y] == '*' else 'X'
        new_board[new_x][new_y] = '*'
        new_board[box_x][box_y] = '$' if board[box_x][box_y] == 'X' else '@'
        
        state = (''.join(''.join(row) for row in new_board))
        if state in visited:
            return None
        visited.add(state)
        return (new_board, (new_x, new_y))
    
    return None

def is_solved(board):
    for row in board:
        for cell in row:
            if cell == '@' or cell == 'X':
                return False
    return True

def solve_sokoban(board, max_moves=30):
    rows, cols = len(board), len(board[0])
    board = [list(row) for row in board]
    
    # Find player position
    player_pos = None
    for i in range(rows):
        for j in range(cols):
            if board[i][j] == '*':
                player_pos = (i, j)
                break
        if player_pos:
            break
    
    def dfs(current_board, pos, moves="", depth=0, visited=None):
        if depth >= max_moves:
            return None
        
        if is_solved(current_board):
            return moves
        
        if visited is None:
            visited = set()
        
        for direction in 'ULDR':  # Try different move orders
            result = try_move(current_board, pos, direction, depth, max_moves, visited)
            if result:
                new_board, new_pos = result
                solution = dfs(new_board, new_pos, moves + direction, depth + 1, visited)
                if solution:
                    return solution
        return None

    return dfs(board, player_pos)

# Initial board
board = [
    "++++++++",
    "+-----X+",
    "+-@X@--+",
    "+X-*X--+",
    "+-@XX@-+",
    "+@--@--+",
    "+------+",
    "++++++++"
]

solution = solve_sokoban(board)
print(solution)