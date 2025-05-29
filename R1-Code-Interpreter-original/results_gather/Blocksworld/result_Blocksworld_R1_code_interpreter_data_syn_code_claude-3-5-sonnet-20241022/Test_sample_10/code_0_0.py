from collections import deque
import copy

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
    new_stacks = copy.deepcopy(stacks)
    block = new_stacks[from_stack-1].pop()
    new_stacks[to_stack-1].append(block)
    return new_stacks

def state_to_string(stacks):
    return str(stacks)

def get_moves(current, target):
    initial_stacks = current
    target_stacks = target
    
    queue = deque([(initial_stacks, [])])
    visited = set()
    
    while queue:
        current_stacks, moves = queue.popleft()
        current_state = state_to_string(current_stacks)
        
        if current_state == state_to_string(target_stacks):
            return moves
            
        if current_state in visited:
            continue
            
        visited.add(current_state)
        
        # Try all possible moves
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack:
                    new_stacks = make_move(current_stacks, from_stack, to_stack)
                    if new_stacks:
                        block = current_stacks[from_stack-1][-1]
                        move = f"Move {block} from {from_stack} to {to_stack}"
                        queue.append((new_stacks, moves + [move]))

# Initial and goal states
initial = [['J'], ['D', 'H', 'F', 'I', 'B', 'A', 'G'], ['E', 'C']]
goal = [['G', 'I'], ['A', 'B', 'C', 'D', 'E', 'F', 'H'], ['J']]

# Get and print solution
solution = get_moves(initial, goal)
print("<<<" + "\n".join(solution) + ">>>")