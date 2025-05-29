from collections import deque
import copy

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return block

def get_state_string(stacks):
    return '|'.join([''.join(stack) for stack in stacks])

def find_solution():
    # Initial state
    initial_stacks = [
        ['G', 'A', 'F', 'E', 'C'],
        ['J', 'H', 'B'],
        ['I', 'D']
    ]
    
    # Goal state
    goal_stacks = [
        ['A', 'B', 'D', 'F', 'I'],
        ['C', 'H'],
        ['E', 'G', 'J']
    ]
    
    # BFS queue: (stacks, moves)
    queue = deque([(initial_stacks, [])])
    visited = {get_state_string(initial_stacks)}
    
    while queue:
        current_stacks, moves = queue.popleft()
        
        # Check if goal reached
        if all(current_stacks[i] == goal_stacks[i] for i in range(3)):
            return moves
        
        # Try all possible moves
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack and is_valid_move(current_stacks, from_stack, to_stack):
                    # Create new state
                    new_stacks = [stack.copy() for stack in current_stacks]
                    block = make_move(new_stacks, from_stack, to_stack)
                    
                    state_string = get_state_string(new_stacks)
                    if state_string not in visited:
                        visited.add(state_string)
                        new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                        queue.append((new_stacks, new_moves))
    
    return None

# Find and print solution
solution = find_solution()
if solution:
    print('\n'.join(solution))
else:
    print("No solution found")