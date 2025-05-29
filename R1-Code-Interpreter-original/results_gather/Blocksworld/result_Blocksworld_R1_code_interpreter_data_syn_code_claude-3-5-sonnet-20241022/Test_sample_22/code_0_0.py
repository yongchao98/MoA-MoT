from collections import deque
import copy

def get_top_blocks(stacks):
    return [stack[-1] if stack else None for stack in stacks]

def is_valid_move(stacks, from_stack, to_stack):
    if from_stack < 0 or from_stack >= len(stacks) or not stacks[from_stack]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
    new_stacks = copy.deepcopy(stacks)
    block = new_stacks[from_stack].pop()
    new_stacks[to_stack].append(block)
    return new_stacks

def state_to_string(stacks):
    return ';'.join(','.join(stack) for stack in stacks)

def find_solution(initial_stacks, goal_stacks):
    queue = deque([(initial_stacks, [])])
    visited = {state_to_string(initial_stacks)}
    
    while queue:
        current_stacks, moves = queue.popleft()
        
        if current_stacks == goal_stacks:
            return moves
        
        for from_stack in range(len(current_stacks)):
            for to_stack in range(len(current_stacks)):
                if from_stack != to_stack:
                    new_stacks = make_move(current_stacks, from_stack, to_stack)
                    if new_stacks:
                        state_str = state_to_string(new_stacks)
                        if state_str not in visited:
                            visited.add(state_str)
                            block = current_stacks[from_stack][-1]
                            new_moves = moves + [f"Move {block} from {from_stack + 1} to {to_stack + 1}"]
                            queue.append((new_stacks, new_moves))
    return None

# Initial and goal states
initial = [['I', 'A', 'G'], ['D', 'F', 'E'], ['H', 'B', 'C']]
goal = [['B', 'C', 'H'], ['F', 'G', 'I'], ['A', 'D', 'E']]

# Find solution
solution = find_solution(initial, goal)

if solution:
    print("<<<")
    print('\n'.join(solution))
    print(">>>")
else:
    print("No solution found")