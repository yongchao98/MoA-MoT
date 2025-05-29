from copy import deepcopy

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
    new_state = deepcopy(state)
    block = new_state[from_stack].pop()
    new_state[to_stack].append(block)
    return new_state

def state_to_string(state):
    return str(state)

def is_goal_state(state, goal):
    return state == goal

def find_solution():
    initial_state = {
        1: ['D', 'A'],
        2: ['B', 'C'],
        3: [],
        4: ['E']
    }
    
    goal_state = {
        1: [],
        2: ['A', 'C', 'D', 'E'],
        3: ['B'],
        4: []
    }
    
    visited = set()
    queue = [(initial_state, [])]
    
    while queue:
        current_state, moves = queue.pop(0)
        state_str = state_to_string(current_state)
        
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        if is_goal_state(current_state, goal_state):
            return moves
            
        for from_stack in range(1, 5):
            if not current_state[from_stack]:
                continue
                
            for to_stack in range(1, 5):
                if from_stack == to_stack:
                    continue
                    
                new_state = make_move(current_state, from_stack, to_stack)
                if new_state:
                    block = current_state[from_stack][-1]
                    new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                    queue.append((new_state, new_moves))

solution = find_solution()
print("<<<" + "\n".join(solution) + ">>>")