from collections import deque

def is_target_state(state, target):
    return target in state

def get_possible_moves(state, capacities):
    moves = []
    for i in range(3):
        # Fill operation
        if state[i] < capacities[i]:
            new_state = list(state)
            new_state[i] = capacities[i]
            moves.append((new_state, f'fill {chr(65+i)}'))
        
        # Empty operation
        if state[i] > 0:
            new_state = list(state)
            new_state[i] = 0
            moves.append((new_state, f'empty {chr(65+i)}'))
        
        # Pour operation
        for j in range(3):
            if i != j and state[i] > 0 and state[j] < capacities[j]:
                new_state = list(state)
                transfer_amount = min(state[i], capacities[j] - state[j])
                new_state[i] -= transfer_amount
                new_state[j] += transfer_amount
                moves.append((new_state, f'pour {chr(65+i)}->{chr(65+j)}'))
    
    return moves

def solve_jug_problem(capacities, target):
    initial_state = (0, 0, 0)
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_target_state(current_state, target):
            return path
        
        for new_state, move in get_possible_moves(current_state, capacities):
            if new_state not in visited:
                visited.add(tuple(new_state))
                queue.append((tuple(new_state), path + [move]))
    
    return None

capacities = (8, 11, 11)
target = 1
solution = solve_jug_problem(capacities, target)
print(solution)