from collections import deque

# Define the initial state of the puzzle
initial_state = {
    'player': (1, 4),
    'boxes': {(2, 2), (2, 3), (3, 5)},
    'goals': {(1, 2), (4, 4)},
    'walls': {(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8),
              (1, 0), (1, 8),
              (2, 0), (2, 5), (2, 8),
              (3, 0), (3, 5), (3, 8),
              (4, 0), (4, 8),
              (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8)}
}

# Define possible moves
moves = {
    'L': (0, -1),
    'R': (0, 1),
    'U': (-1, 0),
    'D': (1, 0)
}

def is_goal_state(state):
    return state['boxes'] == state['goals']

def get_new_position(position, move):
    return (position[0] + move[0], position[1] + move[1])

def is_valid_position(position, state):
    return position not in state['walls'] and position not in state['boxes']

def bfs_solve(initial_state):
    queue = deque([(initial_state, "")])
    visited = set()
    visited.add((initial_state['player'], frozenset(initial_state['boxes'])))

    while queue:
        current_state, path = queue.popleft()

        if is_goal_state(current_state):
            return path

        for move, delta in moves.items():
            new_player_pos = get_new_position(current_state['player'], delta)

            if new_player_pos in current_state['boxes']:
                # If the player is pushing a box
                new_box_pos = get_new_position(new_player_pos, delta)
                if is_valid_position(new_box_pos, current_state):
                    new_boxes = set(current_state['boxes'])
                    new_boxes.remove(new_player_pos)
                    new_boxes.add(new_box_pos)
                    new_state = {
                        'player': new_player_pos,
                        'boxes': frozenset(new_boxes),
                        'goals': current_state['goals'],
                        'walls': current_state['walls']
                    }
                    state_signature = (new_state['player'], new_state['boxes'])
                    if state_signature not in visited:
                        visited.add(state_signature)
                        queue.append((new_state, path + move))
            elif is_valid_position(new_player_pos, current_state):
                # If the player is moving to an empty space
                new_state = {
                    'player': new_player_pos,
                    'boxes': current_state['boxes'],
                    'goals': current_state['goals'],
                    'walls': current_state['walls']
                }
                state_signature = (new_state['player'], new_state['boxes'])
                if state_signature not in visited:
                    visited.add(state_signature)
                    queue.append((new_state, path + move))

    return None

solution = bfs_solve(initial_state)
print(solution)