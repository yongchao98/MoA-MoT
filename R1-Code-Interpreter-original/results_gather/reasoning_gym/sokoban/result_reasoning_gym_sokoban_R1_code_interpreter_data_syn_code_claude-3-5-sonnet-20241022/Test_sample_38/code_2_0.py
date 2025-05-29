from collections import deque

def is_valid(x, y, rows, cols):
    return 0 <= x < rows and 0 <= y < cols

def get_next_states(state, rows, cols):
    moves = []
    player1, player2 = None, None
    boxes = []
    goals = []
    
    # Find players, boxes and goals
    for i in range(rows):
        for j in range(cols):
            if state[i][j] == '*':
                player1 = (i, j)
            elif state[i][j] == '@':
                if not player2:  # First @ encountered is player2
                    player2 = (i, j)
                else:  # Other @ are boxes
                    boxes.append((i, j))
            elif state[i][j] == '$':
                boxes.append((i, j))
            elif state[i][j] == 'X':
                goals.append((i, j))

    directions = [(0, 1, 'R'), (0, -1, 'L'), (1, 0, 'D'), (-1, 0, 'U')]
    
    # Try moves for both players
    for player in [player1, player2]:
        px, py = player
        for dx, dy, move in directions:
            new_x, new_y = px + dx, py + dy
            
            if not is_valid(new_x, new_y, rows, cols):
                continue
                
            new_state = [list(row) for row in state]
            curr_pos = new_state[px][py]
            
            # Moving to empty space
            if new_state[new_x][new_y] in ['-', 'X']:
                new_state[new_x][new_y] = '*'
                new_state[px][py] = '-'
                moves.append((''.join(''.join(row) for row in new_state), move))
                
            # Moving a box
            elif new_state[new_x][new_y] in ['@', '$']:
                push_x, push_y = new_x + dx, new_y + dy
                if (is_valid(push_x, push_y, rows, cols) and 
                    new_state[push_x][push_y] in ['-', 'X']):
                    # Update positions
                    new_state[push_x][push_y] = '@' if new_state[push_x][push_y] == '-' else '$'
                    new_state[new_x][new_y] = '*'
                    new_state[px][py] = '-'
                    moves.append((''.join(''.join(row) for row in new_state), move))
    
    return moves

def is_goal(state):
    goals = 0
    boxes_on_goals = 0
    for row in state:
        for cell in row:
            if cell in ['X', '$']:
                goals += 1
            if cell == '$':
                boxes_on_goals += 1
    return goals == 3 and boxes_on_goals == 3

def solve_sokoban(initial_state):
    initial = [''.join(row) for row in initial_state]
    initial = ''.join(initial)
    
    queue = deque([(initial, "")])
    visited = {initial}
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal(current_state):
            return path
            
        current_state_2d = [list(current_state[i:i+6]) for i in range(0, len(current_state), 6)]
        
        for next_state, move in get_next_states(current_state_2d, 8, 6):
            if next_state not in visited:
                visited.add(next_state)
                queue.append((next_state, path + move))
    
    return None

initial_state = [
    ['+', '+', '+', '+', '+', '+'],
    ['+', '+', 'X', '@', '*', '+'],
    ['+', '+', '+', '$', '@', '+'],
    ['+', '+', '+', '+', '-', '+'],
    ['+', '+', '$', '+', '-', '+'],
    ['+', '+', '+', '+', '-', '+'],
    ['+', '+', '$', '+', 'X', '+'],
    ['+', '+', '+', '+', '+', '+']
]

solution = solve_sokoban(initial_state)
print(solution)