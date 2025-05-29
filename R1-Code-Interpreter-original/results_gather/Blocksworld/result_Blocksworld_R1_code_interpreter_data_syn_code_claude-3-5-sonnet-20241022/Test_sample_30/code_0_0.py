from collections import deque
import copy

def get_top_blocks(stacks):
    top_blocks = []
    for i, stack in enumerate(stacks):
        if stack:
            top_blocks.append((stack[-1], i))
    return top_blocks

def is_goal_state(current, goal):
    return all(current[i] == goal[i] for i in range(len(current)))

def state_to_string(state):
    return ';'.join(','.join(stack) for stack in state)

def find_solution():
    initial = [['C', 'F'], ['D'], ['A', 'B', 'E']]
    goal = [['E', 'F'], ['A', 'B', 'C', 'D'], []]
    
    queue = deque([(initial, [])])
    visited = {state_to_string(initial)}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if is_goal_state(current_state, goal):
            return moves
        
        # Get all possible moves from current state
        top_blocks = get_top_blocks(current_state)
        
        for block, source in top_blocks:
            for dest in range(3):
                if dest != source:
                    # Try moving block from source to dest
                    new_state = copy.deepcopy(current_state)
                    new_state[source].pop()
                    new_state[dest].append(block)
                    
                    state_str = state_to_string(new_state)
                    if state_str not in visited:
                        visited.add(state_str)
                        move = f"Move {block} from {source+1} to {dest+1}"
                        queue.append((new_state, moves + [move]))

solution = find_solution()
print('\n'.join(solution))