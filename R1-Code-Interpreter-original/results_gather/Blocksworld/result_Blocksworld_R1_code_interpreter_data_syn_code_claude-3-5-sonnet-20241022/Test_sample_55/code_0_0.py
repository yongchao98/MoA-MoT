from collections import deque
import copy

def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack has blocks
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return stacks

def state_to_string(stacks):
    return ';'.join(','.join(stack) for stack in stacks)

def is_goal_state(current, goal):
    return all(current[i] == goal[i] for i in range(len(current)))

def find_solution():
    # Initial state
    initial_stacks = [
        ['J'],
        ['C', 'K', 'F', 'E', 'B', 'I'],
        ['G', 'A'],
        ['H', 'D']
    ]
    
    # Goal state
    goal_stacks = [
        ['F', 'G', 'H', 'I', 'J'],
        ['C', 'D', 'E', 'K'],
        ['A', 'B'],
        []
    ]
    
    queue = deque([(initial_stacks, [])])
    visited = {state_to_string(initial_stacks)}
    
    while queue:
        current_stacks, moves = queue.popleft()
        
        if is_goal_state(current_stacks, goal_stacks):
            return moves
        
        # Try all possible moves
        for from_stack in range(1, 5):
            for to_stack in range(1, 5):
                if from_stack != to_stack:
                    if is_valid_move(current_stacks, from_stack, to_stack):
                        # Create a copy of current state
                        new_stacks = [stack[:] for stack in current_stacks]
                        block = get_top_block(new_stacks[from_stack-1])
                        make_move(new_stacks, from_stack, to_stack)
                        
                        new_state = state_to_string(new_stacks)
                        if new_state not in visited:
                            visited.add(new_state)
                            new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                            queue.append((new_stacks, new_moves))
    
    return None

# Find and print solution
solution = find_solution()
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")