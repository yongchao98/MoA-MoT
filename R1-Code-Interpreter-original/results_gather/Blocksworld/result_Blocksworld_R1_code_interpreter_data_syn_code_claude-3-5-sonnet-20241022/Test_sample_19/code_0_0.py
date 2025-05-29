from collections import deque
import copy

def get_top_blocks(state):
    return [stack[-1] if stack else None for stack in state]

def is_valid_move(state, from_stack, to_stack):
    # Check if source stack has blocks and destination exists
    if not state[from_stack] or to_stack >= len(state):
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
    return state_to_string(current) == state_to_string(goal)

def find_solution():
    initial_state = [['G', 'F'], ['E'], ['A', 'B'], ['D', 'C']]
    goal_state = [['B', 'D', 'E'], ['G'], ['A', 'C', 'F']]
    
    queue = deque([(initial_state, [])])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if is_goal_state(current_state, goal_state):
            return moves
        
        for from_stack in range(len(current_state)):
            for to_stack in range(len(current_state)):
                if from_stack != to_stack:
                    if is_valid_move(current_state, from_stack, to_stack):
                        new_state = make_move(current_state, from_stack, to_stack)
                        state_str = state_to_string(new_state)
                        
                        if state_str not in visited:
                            visited.add(state_str)
                            block = current_state[from_stack][-1]
                            move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                            queue.append((new_state, moves + [move]))
    
    return None

solution = find_solution()
if solution:
    print("<<<")
    print('\n'.join(solution))
    print(">>>")
else:
    print("No solution found")