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
    return block

def get_state_string(stacks):
    return str(stacks)

def is_goal_state(current, goal):
    return current == goal

def find_solution():
    # Initial state
    initial_stacks = [['E', 'A', 'B', 'C'], [], ['D', 'F']]
    # Goal state
    goal_stacks = [['A', 'C', 'E'], ['B', 'F'], ['D']]
    
    # BFS queue: (stacks, moves)
    queue = deque([(initial_stacks, [])])
    visited = {get_state_string(initial_stacks)}
    
    while queue:
        current_stacks, moves = queue.popleft()
        
        if is_goal_state(current_stacks, goal_stacks):
            return moves
        
        # Try all possible moves
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack:
                    # Create a copy of current state
                    new_stacks = [stack[:] for stack in current_stacks]
                    
                    if is_valid_move(new_stacks, from_stack, to_stack):
                        block = make_move(new_stacks, from_stack, to_stack)
                        new_state = get_state_string(new_stacks)
                        
                        if new_state not in visited:
                            visited.add(new_state)
                            new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                            queue.append((new_stacks, new_moves))
    
    return None

# Find and print solution
solution = find_solution()
if solution:
    print("<<<")
    print("\n".join(solution))
    print(">>>")
else:
    print("No solution found")