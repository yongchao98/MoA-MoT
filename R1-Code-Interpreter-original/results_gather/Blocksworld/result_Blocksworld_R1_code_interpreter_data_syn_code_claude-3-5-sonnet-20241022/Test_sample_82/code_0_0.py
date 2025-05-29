from collections import deque
import copy

def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(stacks, from_stack, to_stack):
    if not stacks[from_stack-1]:  # Empty source stack
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
    return ';'.join(','.join(stack) for stack in stacks)

def is_goal_state(current, goal):
    return all(current[i] == goal[i] for i in range(len(current)))

def find_solution():
    initial = [['F', 'A'], ['K', 'D', 'E', 'B', 'I', 'C'], ['H', 'J', 'G']]
    goal = [['G'], ['B', 'D'], ['A', 'C', 'E', 'F', 'H', 'I', 'J', 'K']]
    
    queue = deque([(initial, [])])
    visited = {state_to_string(initial)}
    
    while queue:
        current_stacks, moves = queue.popleft()
        
        if is_goal_state(current_stacks, goal):
            return moves
        
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack:
                    new_stacks = make_move(current_stacks, from_stack, to_stack)
                    if new_stacks:
                        state_str = state_to_string(new_stacks)
                        if state_str not in visited:
                            visited.add(state_str)
                            block = current_stacks[from_stack-1][-1]
                            new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                            queue.append((new_stacks, new_moves))
    return None

solution = find_solution()
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")