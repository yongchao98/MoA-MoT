from collections import deque
import copy

def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack has blocks
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, move_from, move_to):
    if not is_valid_move(stacks, move_from, move_to):
        return None
    new_stacks = copy.deepcopy(stacks)
    block = new_stacks[move_from-1].pop()
    new_stacks[move_to-1].append(block)
    return new_stacks

def state_to_string(stacks):
    return ';'.join(','.join(stack) for stack in stacks)

def is_goal_state(current, goal):
    return all(current[i] == goal[i] for i in range(len(current)))

def find_solution():
    initial = [['F', 'B', 'I', 'K', 'A'], ['H', 'D'], ['J', 'E', 'C', 'G']]
    goal = [['C', 'G'], ['A', 'B', 'D', 'E', 'I', 'K'], ['F', 'H', 'J']]
    
    queue = deque([(initial, [])])
    seen = {state_to_string(initial)}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if is_goal_state(current_state, goal):
            return moves
        
        # Try all possible moves
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack:
                    new_state = make_move(current_state, from_stack, to_stack)
                    if new_state:
                        state_str = state_to_string(new_state)
                        if state_str not in seen:
                            block = current_state[from_stack-1][-1]
                            moves_copy = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                            queue.append((new_state, moves_copy))
                            seen.add(state_str)
    return None

solution = find_solution()
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")