from collections import deque
import copy

def get_top_blocks(state):
    return [stack[-1] if stack else None for stack in state]

def is_valid_move(state, from_stack, to_stack):
    # Check if source stack has blocks and destination stack has space
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

def find_solution(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if current_state == goal_state:
            return moves
            
        for from_stack in range(len(current_state)):
            for to_stack in range(len(current_state)):
                if from_stack != to_stack and current_state[from_stack]:
                    new_state = make_move(current_state, from_stack, to_stack)
                    if new_state:
                        state_str = state_to_string(new_state)
                        if state_str not in visited:
                            visited.add(state_str)
                            move = f"Move {current_state[from_stack][-1]} from {from_stack + 1} to {to_stack + 1}"
                            queue.append((new_state, moves + [move]))
    return None

# Initial and goal states
initial_state = [
    ['C', 'E'],
    ['D'],
    ['G', 'F', 'B'],
    ['A', 'H']
]

goal_state = [
    ['B', 'C', 'D', 'E', 'F', 'G', 'H'],
    [],
    ['A'],
    []
]

# Find solution
solution = find_solution(initial_state, goal_state)

if solution:
    print("<<<")
    print('\n'.join(solution))
    print(">>>")
else:
    print("No solution found")