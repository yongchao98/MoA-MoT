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
    # Get the block from top of source stack
    block = stacks[from_stack-1].pop()
    # Add it to destination stack
    stacks[to_stack-1].append(block)
    return stacks

def get_state_string(stacks):
    return str(stacks)

def find_solution(initial_stacks, goal_stacks):
    visited = set()
    queue = deque([(initial_stacks, [])])
    
    while queue:
        current_stacks, moves = queue.popleft()
        current_state = get_state_string(current_stacks)
        
        if current_state == get_state_string(goal_stacks):
            return moves
            
        if current_state in visited:
            continue
            
        visited.add(current_state)
        
        # Try all possible moves
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack:
                    if not current_stacks[from_stack-1]:  # Skip if source stack is empty
                        continue
                    
                    new_stacks = copy.deepcopy(current_stacks)
                    block = new_stacks[from_stack-1][-1]  # Get top block
                    new_stacks = make_move(new_stacks, from_stack, to_stack)
                    
                    if new_stacks:
                        new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                        queue.append((new_stacks, new_moves))

# Initial and goal states
initial = [[], ['B'], ['C', 'A', 'E', 'D']]
goal = [['D'], ['C'], ['A', 'B', 'E']]

# Find solution
solution = find_solution(initial, goal)

# Print solution
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")