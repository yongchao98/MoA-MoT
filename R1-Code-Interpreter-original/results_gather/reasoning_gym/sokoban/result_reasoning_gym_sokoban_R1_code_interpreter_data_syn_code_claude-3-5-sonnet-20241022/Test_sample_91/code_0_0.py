from collections import deque
import copy

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_state_string(board):
    return '\n'.join([''.join(row) for row in board])

def is_win(board):
    # Check if all boxes are on goals
    for row in board:
        for cell in row:
            if cell == '@':  # If there's still a box not on a goal
                return False
    return True

def solve_sokoban(initial_board):
    rows = len(initial_board)
    cols = len(initial_board[0])
    
    # Convert initial board to list of lists for mutability
    board = [list(row) for row in initial_board]
    
    # Directions: Up, Down, Left, Right
    directions = [(-1, 0, 'U'), (1, 0, 'D'), (0, -1, 'L'), (0, 1, 'R')]
    
    visited = set()
    queue = deque([(board, "", [])])  # (board, moves, player_positions)
    
    while queue:
        current_board, moves, players = queue.popleft()
        
        if not players:  # First time, find all players
            players = []
            for i in range(rows):
                for j in range(cols):
                    if current_board[i][j] in ['@', '*']:
                        players.append((i, j))
        
        current_state = get_state_string(current_board)
        if current_state in visited:
            continue
            
        visited.add(current_state)
        
        if is_win(current_board):
            return moves
            
        # Try moving each player
        for player_pos in players:
            px, py = player_pos
            
            for dx, dy, move in directions:
                new_x, new_y = px + dx, py + dy
                
                if not is_valid(new_x, new_y, rows, cols):
                    continue
                    
                # If moving to empty space or goal
                if current_board[new_x][new_y] in ['-', 'X']:
                    new_board = [row[:] for row in current_board]
                    new_players = players.copy()
                    
                    # Update player position
                    new_players.remove((px, py))
                    new_players.append((new_x, new_y))
                    
                    # Update board
                    if current_board[px][py] == '*':
                        new_board[px][py] = 'X'
                    else:
                        new_board[px][py] = '-'
                        
                    if current_board[new_x][new_y] == 'X':
                        new_board[new_x][new_y] = '*'
                    else:
                        new_board[new_x][new_y] = '@'
                    
                    queue.append((new_board, moves + move, new_players))

# Initial board
initial = [
    "++++++",
    "+-@-X+",
    "+-$-++",
    "+X@X-+",
    "+-@--+",
    "+-@X*+",
    "+----+",
    "++++++",
]

solution = solve_sokoban(initial)
print(solution)