from collections import deque
import copy

def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(state, from_stack, to_stack):
    # Check if source stack has blocks
    if not state[from_stack-1]:
        return False
    return True

def make_move(state, move_from, move_to):
    if not is_valid_move(state, move_from, move_to):
        return None
    new_state = copy.deepcopy(state)
    block = new_state[move_from-1].pop()
    new_state[move_to-1].append(block)
    return new_state

def state_to_string(state):
    return ','.join([''.join(stack) for stack in state])

def find_solution(initial_state, goal_state):
    visited = set()
    queue = deque([(initial_state, [])])
    visited.add(state_to_string(initial_state))
    
    while queue:
        current_state, moves = queue.popleft()
        
        if current_state == goal_state:
            return moves
            
        # Try all possible moves
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack and current_state[from_stack-1]:
                    new_state = make_move(current_state, from_stack, to_stack)
                    if new_state:
                        state_str = state_to_string(new_state)
                        if state_str not in visited:
                            block = current_state[from_stack-1][-1]
                            move = f"Move {block} from {from_stack} to {to_stack}"
                            queue.append((new_state, moves + [move]))
                            visited.add(state_str)
    return None

# Initial and goal states
initial = [['G', 'C'], ['D', 'A', 'E', 'B'], ['I', 'H', 'F']]
goal = [['B', 'I'], ['E'], ['A', 'C', 'D', 'F', 'G', 'H']]

# Find solution
solution = find_solution(initial, goal)
if solution:
    print('<<<')
    print('\n'.join(solution))
    print('>>>')
else:
    print("No solution found")