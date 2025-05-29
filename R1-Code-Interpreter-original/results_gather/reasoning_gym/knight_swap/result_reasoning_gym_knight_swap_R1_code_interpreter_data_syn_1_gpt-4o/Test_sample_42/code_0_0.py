from collections import deque

# Define the board positions
positions = {
    'A3': (0, 0), 'B3': (0, 1), 'C3': (0, 2), 'D3': (0, 3),
    'A2': (1, 0), 'B2': (1, 1), 'C2': (1, 2), 'D2': (1, 3),
    'A1': (2, 0), 'B1': (2, 1), 'C1': (2, 2), 'D1': (2, 3)
}

# Reverse lookup for positions
reverse_positions = {v: k for k, v in positions.items()}

# Initial and target positions
initial_state = {'w': [(0, 1), (0, 2)], 'B': [(0, 0), (1, 3)]}
target_state = {'w': [(0, 0), (1, 3)], 'B': [(0, 1), (0, 2)]}

# Knight moves (L-shape)
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1),
                (1, 2), (1, -2), (-1, 2), (-1, -2)]

def is_valid_position(pos):
    return 0 <= pos[0] < 3 and 0 <= pos[1] < 4

def swap_knights():
    queue = deque([(initial_state, [], 'w')])  # (current_state, moves, current_turn)
    visited = set()

    while queue:
        current_state, moves, turn = queue.popleft()
        state_key = (tuple(sorted(current_state['w'])), tuple(sorted(current_state['B'])))
        
        if state_key in visited:
            continue
        visited.add(state_key)

        # Check if we reached the target state
        if current_state == target_state:
            return moves

        # Get the positions of the current player's knights
        for i, knight_pos in enumerate(current_state[turn]):
            for move in knight_moves:
                new_pos = (knight_pos[0] + move[0], knight_pos[1] + move[1])
                if is_valid_position(new_pos) and new_pos not in current_state['w'] and new_pos not in current_state['B']:
                    # Create a new state with the knight moved
                    new_state = {k: v[:] for k, v in current_state.items()}
                    new_state[turn][i] = new_pos
                    new_moves = moves + [f"{turn},{reverse_positions[knight_pos]},{reverse_positions[new_pos]}"]
                    # Alternate turn
                    next_turn = 'B' if turn == 'w' else 'w'
                    queue.append((new_state, new_moves, next_turn))

    return "No"

result = swap_knights()
print(result)