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
    return "|".join(",".join(stack) for stack in state)

def is_goal_state(state):
    goal = [
        ['A', 'B', 'C', 'D', 'G', 'H', 'I'],
        ['E'],
        ['F'],
        []
    ]
    return state == goal

def solve_blocksworld():
    initial_state = [
        ['G', 'A', 'H'],
        ['E'],
        ['I', 'C', 'B'],
        ['F', 'D']
    ]
    
    visited = set()
    queue = deque([(initial_state, [])])
    
    while queue:
        current_state, moves = queue.popleft()
        state_str = state_to_string(current_state)
        
        if is_goal_state(current_state):
            return moves
        
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        # Try all possible moves
        for from_stack in range(4):
            if not current_state[from_stack]:
                continue
            
            block = current_state[from_stack][-1]
            for to_stack in range(4):
                if from_stack != to_stack:
                    new_state = make_move(current_state, from_stack, to_stack)
                    if new_state:
                        move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                        queue.append((new_state, moves + [move]))

solution = solve_blocksworld()
print("<<<" + "\n".join(solution) + ">>>")