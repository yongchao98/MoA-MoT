from collections import deque
import copy

def is_valid_move(state, from_stack, to_stack):
    # Check if source stack is not empty
    if not state[from_stack]:
        return False
    return True

def make_move(state, from_stack, to_stack):
    if not is_valid_move(state, from_stack, to_stack):
        return None
    new_state = copy.deepcopy(state)
    block = new_state[from_stack].pop()
    new_state[to_stack].append(block)
    return new_state

def state_to_string(state):
    return str(state)

def is_goal_state(state, goal):
    return state == goal

def find_solution():
    # Initial state
    initial_state = {
        1: ['D', 'B'],
        2: [],
        3: ['C'],
        4: ['A', 'E', 'F']
    }
    
    # Goal state
    goal_state = {
        1: ['B', 'D'],
        2: ['A', 'C', 'E', 'F'],
        3: []
    }
    
    queue = deque([(initial_state, [])])
    visited = set()
    
    while queue:
        current_state, path = queue.popleft()
        state_str = state_to_string(current_state)
        
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        if is_goal_state(current_state, goal_state):
            return path
        
        # Try all possible moves
        for from_stack in range(1, 5):
            for to_stack in range(1, 5):
                if from_stack != to_stack and current_state[from_stack]:
                    new_state = make_move(current_state, from_stack, to_stack)
                    if new_state:
                        block = current_state[from_stack][-1]
                        move = f"Move {block} from {from_stack} to {to_stack}"
                        queue.append((new_state, path + [move]))
    
    return None

solution = find_solution()
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")