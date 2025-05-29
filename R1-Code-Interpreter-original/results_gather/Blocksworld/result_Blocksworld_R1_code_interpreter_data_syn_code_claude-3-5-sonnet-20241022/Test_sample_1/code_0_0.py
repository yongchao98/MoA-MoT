from collections import deque
import copy

def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(state, from_stack, to_stack):
    # Check if source stack has blocks
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

def find_solution(initial_state, goal_state):
    visited = set()
    queue = deque([(initial_state, [])])
    
    while queue:
        current_state, moves = queue.popleft()
        state_str = state_to_string(current_state)
        
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        if is_goal_state(current_state, goal_state):
            return moves
            
        # Try all possible moves
        for from_stack in range(4):
            for to_stack in range(4):
                if from_stack != to_stack and current_state[from_stack]:
                    new_state = make_move(current_state, from_stack, to_stack)
                    if new_state:
                        block = current_state[from_stack][-1]
                        move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                        queue.append((new_state, moves + [move]))
    
    return None

# Initial state
initial = [
    [],         # stack1
    ['A', 'B'], # stack2
    ['D'],      # stack3
    ['C', 'E']  # stack4
]

# Goal state
goal = [
    [],             # stack1
    ['A', 'D', 'E'], # stack2
    ['B', 'C'],    # stack3
    []              # stack4
]

solution = find_solution(initial, goal)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")