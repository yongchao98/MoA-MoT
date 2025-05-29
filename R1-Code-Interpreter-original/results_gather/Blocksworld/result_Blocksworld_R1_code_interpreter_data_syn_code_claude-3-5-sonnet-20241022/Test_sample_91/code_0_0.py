from collections import deque
import copy

def get_top_blocks(state):
    tops = []
    for i, stack in enumerate(state):
        if stack:
            tops.append((stack[-1], i+1))
    return tops

def is_valid_move(from_stack, to_stack, state):
    # Check if from_stack has blocks
    if not state[from_stack-1]:
        return False
    return True

def make_move(state, from_stack, to_stack):
    new_state = copy.deepcopy(state)
    block = new_state[from_stack-1].pop()
    new_state[to_stack-1].append(block)
    return new_state

def state_to_string(state):
    return str(state)

def is_goal_state(state, goal):
    return state == goal

def find_solution():
    initial = [['E'], ['B'], ['H', 'C', 'A', 'F'], ['D', 'G']]
    goal = [['D', 'E', 'G', 'H'], ['B', 'C', 'F'], ['A'], []]
    
    visited = set()
    queue = deque([(initial, [])])
    
    while queue:
        current_state, moves = queue.popleft()
        state_str = state_to_string(current_state)
        
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        if is_goal_state(current_state, goal):
            return moves
            
        # Get all possible moves
        top_blocks = get_top_blocks(current_state)
        
        for block, from_stack in top_blocks:
            for to_stack in range(1, 5):  # 4 stacks
                if from_stack != to_stack and is_valid_move(from_stack, to_stack, current_state):
                    new_state = make_move(current_state, from_stack, to_stack)
                    new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                    queue.append((new_state, new_moves))

solution = find_solution()
print("<<<" + "\n".join(solution) + ">>>")