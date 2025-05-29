from collections import deque
import copy

def is_valid(x, y, width, height):
    return 0 <= x < height and 0 <= y < width

def get_next_states(state):
    directions = [('U', -1, 0), ('D', 1, 0), ('L', 0, -1), ('R', 0, 1)]
    next_states = []
    height = len(state['board'])
    width = len(state['board'][0])
    px, py = state['player']
    
    for direction, dx, dy in directions:
        new_x, new_y = px + dx, py + dy
        
        if not is_valid(new_x, new_y, width, height) or state['board'][new_x][new_y] == '+':
            continue
            
        new_state = {'board': [row[:] for row in state['board']], 'player': (new_x, new_y), 'path': state['path'] + direction}
        
        # If moving to empty space or goal
        if state['board'][new_x][new_y] in '-X':
            new_state['board'][px][py] = '-' if state['board'][px][py] in '*%' else 'X'
            new_state['board'][new_x][new_y] = '*' if state['board'][new_x][new_y] == '-' else '%'
            next_states.append(new_state)
            
        # If moving a box
        elif state['board'][new_x][new_y] in '@$':
            push_x, push_y = new_x + dx, new_y + dy
            if is_valid(push_x, push_y, width, height) and state['board'][push_x][push_y] in '-X':
                new_state['board'][px][py] = '-' if state['board'][px][py] in '*%' else 'X'
                new_state['board'][new_x][new_y] = '*' if state['board'][new_x][new_y] == '@' else '%'
                new_state['board'][push_x][push_y] = '@' if state['board'][push_x][push_y] == '-' else '$'
                next_states.append(new_state)
                
    return next_states

def is_solved(board):
    for row in board:
        for cell in row:
            if cell == '@' or cell == 'X':
                return False
    return True

def solve_sokoban():
    initial_board = [
        ['+', '+', '+', '+', '+', '+'],
        ['+', '-', '-', '-', '-', '+'],
        ['+', '-', 'X', '-', '-', '+'],
        ['+', '@', '@', '-', '-', '+'],
        ['+', '%', '@', 'X', '$', '+'],
        ['+', '+', '+', '+', '+', '+']
    ]
    
    initial_state = {
        'board': initial_board,
        'player': (4, 1),
        'path': ''
    }
    
    visited = set()
    queue = deque([initial_state])
    
    while queue:
        current_state = queue.popleft()
        board_tuple = tuple(tuple(row) for row in current_state['board'])
        player_pos = current_state['player']
        state_key = (board_tuple, player_pos)
        
        if state_key in visited:
            continue
            
        visited.add(state_key)
        
        if is_solved(current_state['board']):
            return current_state['path']
            
        next_states = get_next_states(current_state)
        queue.extend(next_states)
    
    return None

solution = solve_sokoban()
print(solution)