from collections import deque

def is_valid_state(state):
    goal = [['E', 'F'], ['A', 'B', 'C', 'D'], []]
    # Check if blocks in stack1 are in correct order
    if len(state[0]) > 0:
        if len(state[0]) > 2: return False
        if len(state[0]) == 2 and state[0] != ['E', 'F']: return False
        if len(state[0]) == 1 and state[0][0] not in ['E', 'F']: return False
    
    # Check if blocks in stack2 are in correct order
    if len(state[1]) > 0:
        correct_order = ['A', 'B', 'C', 'D']
        for i, block in enumerate(state[1]):
            if i >= len(correct_order): return False
            if block != correct_order[i]: return False
    
    return True

def get_next_states(state):
    moves = []
    for src in range(3):
        if not state[src]:  # Skip empty stacks
            continue
        for dst in range(3):
            if src != dst:
                new_state = [stack[:] for stack in state]
                block = new_state[src].pop()
                new_state[dst].append(block)
                if is_valid_state(new_state):
                    moves.append((f"Move {block} from {src+1} to {dst+1}", new_state))
    return moves

def is_goal_state(state):
    return state == [['E', 'F'], ['A', 'B', 'C', 'D'], []]

def solve_blocks():
    initial = [['C', 'F'], ['D'], ['A', 'B', 'E']]
    queue = deque([(initial, [])])
    visited = set()
    
    while queue:
        current_state, moves = queue.popleft()
        state_str = str(current_state)
        
        if state_str in visited:
            continue
        visited.add(state_str)
        
        if is_goal_state(current_state):
            return moves
        
        for move, new_state in get_next_states(current_state):
            if str(new_state) not in visited:
                queue.append((new_state, moves + [move]))
    
    return None

solution = solve_blocks()
if solution:
    print('\n'.join(solution))