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
    return ';'.join(','.join(stack) for stack in state)

def is_goal_state(current, goal):
    return all(current[i] == goal[i] for i in range(len(current)))

def find_solution():
    initial_state = [['J', 'B'], ['F', 'I', 'E'], ['D', 'A', 'C', 'H', 'G']]
    goal_state = [['B', 'C', 'G', 'H', 'I', 'J'], ['D', 'F'], ['A', 'E']]
    
    queue = deque([(initial_state, [])])
    visited = set()
    
    while queue:
        current_state, moves = queue.popleft()
        current_str = state_to_string(current_state)
        
        if current_str in visited:
            continue
        
        visited.add(current_str)
        
        if is_goal_state(current_state, goal_state):
            return moves
        
        # Try all possible moves
        for from_stack in range(3):
            for to_stack in range(3):
                if from_stack != to_stack and current_state[from_stack]:
                    new_state = make_move(current_state, from_stack, to_stack)
                    if new_state:
                        block = current_state[from_stack][-1]
                        move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                        queue.append((new_state, moves + [move]))

solution = find_solution()
print("<<<" + "\n".join(solution) + ">>>")