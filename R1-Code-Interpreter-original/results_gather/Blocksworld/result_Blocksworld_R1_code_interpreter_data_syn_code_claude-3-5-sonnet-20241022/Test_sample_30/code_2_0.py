from collections import deque
import copy

def get_legal_moves(state):
    moves = []
    # For each source stack
    for src in range(3):
        if not state[src]:  # Skip empty stacks
            continue
        # For each destination stack
        for dst in range(3):
            if src != dst:  # Can't move to same stack
                moves.append((src, dst, state[src][-1]))
    return moves

def apply_move(state, src, dst):
    new_state = copy.deepcopy(state)
    block = new_state[src].pop()
    new_state[dst].append(block)
    return new_state

def is_goal_state(state):
    goal = [['E', 'F'], ['A', 'B', 'C', 'D'], []]
    return state == goal

def state_to_string(state):
    return ';'.join(','.join(stack) for stack in state)

def find_solution():
    initial = [['C', 'F'], ['D'], ['A', 'B', 'E']]
    visited = {state_to_string(initial)}
    queue = deque([(initial, [])])
    
    while queue:
        current_state, path = queue.popleft()
        
        if is_goal_state(current_state):
            return path
            
        for src, dst, block in get_legal_moves(current_state):
            new_state = apply_move(current_state, src, dst)
            state_str = state_to_string(new_state)
            
            if state_str not in visited:
                visited.add(state_str)
                move = f"Move {block} from {src+1} to {dst+1}"
                new_path = path + [move]
                queue.append((new_state, new_path))

# Find and print solution
solution = find_solution()
if solution:
    print('\n'.join(solution))