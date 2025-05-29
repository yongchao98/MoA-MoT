from collections import deque

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
            
        new_board = [row[:] for row in state['board']]
        new_path = state['path'] + direction
        
        # Check if current position is player-on-goal
        is_player_on_goal = state['board'][px][py] == '%'
        
        # Moving to empty space or goal
        if state['board'][new_x][new_y] in '-X':
            # Update current position
            new_board[px][py] = 'X' if is_player_on_goal else '-'
            # Update new position
            new_board[new_x][new_y] = '%' if state['board'][new_x][new_y] == 'X' else '*'
            next_states.append({'board': new_board, 'player': (new_x, new_y), 'path': new_path})
            
        # Moving a box
        elif state['board'][new_x][new_y] in '@$':
            push_x, push_y = new_x + dx, new_y + dy
            if is_valid(push_x, push_y, width, height) and state['board'][push_x][push_y] in '-X':
                # Update current position
                new_board[px][py] = 'X' if is_player_on_goal else '-'
                # Update box position
                new_board[new_x][new_y] = '%' if state['board'][new_x][new_y] == '$' else '*'
                # Update pushed box position
                new_board[push_x][push_y] = '$' if state['board'][push_x][push_y] == 'X' else '@'
                next_states.append({'board': new_board, 'player': (new_x, new_y), 'path': new_path})
    
    return next_states

def is_solved(board):
    for row in board:
        for cell in row:
            if cell == '@' or cell == 'X':  # If there's any box not on goal or empty goal
                return False
    return True

def board_to_string(board):
    return ''.join(''.join(row) for row in board)

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
        board_string = board_to_string(current_state['board'])
        player_pos = current_state['player']
        state_key = (board_string, player_pos)
        
        if state_key in visited:
            continue
            
        visited.add(state_key)
        
        if is_solved(current_state['board']):
            return current_state['path']
            
        for next_state in get_next_states(current_state):
            queue.append(next_state)
    
    return None

solution = solve_sokoban()
print(solution)