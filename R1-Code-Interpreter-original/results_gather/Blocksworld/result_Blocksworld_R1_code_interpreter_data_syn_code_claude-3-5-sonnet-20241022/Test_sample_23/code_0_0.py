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

def solve_blocksworld():
    # Initial state
    initial_stacks = [
        ['B', 'I'],
        ['G', 'F'],
        ['C', 'D', 'E'],
        ['A', 'H']
    ]
    
    # Goal state
    goal_stacks = [
        ['E'],
        ['A', 'D'],
        ['B', 'C', 'F', 'G', 'H', 'I']
    ]
    
    goal_state = get_state_string(goal_stacks)
    
    # BFS queue: (stacks, moves)
    queue = deque([(initial_stacks, [])])
    visited = set()
    
    while queue:
        current_stacks, moves = queue.popleft()
        current_state = get_state_string(current_stacks)
        
        if current_state == goal_state:
            return moves
            
        if current_state in visited:
            continue
            
        visited.add(current_state)
        
        # Try all possible moves
        for from_stack in range(1, 5):
            for to_stack in range(1, 5):
                if from_stack != to_stack:
                    if is_valid_move(current_stacks, from_stack, to_stack):
                        # Create a copy of current state
                        new_stacks = [stack.copy() for stack in current_stacks]
                        block = make_move(new_stacks, from_stack, to_stack)
                        new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                        queue.append((new_stacks, new_moves))
    
    return None

# Solve and print the solution
solution = solve_blocksworld()
if solution:
    print('\n'.join(solution))
else:
    print("No solution found")