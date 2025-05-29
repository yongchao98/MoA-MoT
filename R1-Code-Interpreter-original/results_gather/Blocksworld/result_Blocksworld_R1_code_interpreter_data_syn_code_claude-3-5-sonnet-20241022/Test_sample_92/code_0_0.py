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

def get_state_string(stacks):
    return '|'.join([''.join(stack) for stack in stacks])

def find_solution(initial_stacks, goal_stacks):
    queue = deque([(initial_stacks, [])])
    visited = set()
    
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
                    if is_valid_move(current_stacks, from_stack, to_stack):
                        new_stacks = [stack.copy() for stack in current_stacks]
                        block = get_top_block(new_stacks[from_stack-1])
                        make_move(new_stacks, from_stack, to_stack)
                        new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                        queue.append((new_stacks, new_moves))

# Initial and goal states
initial = [['F'], ['G', 'B', 'D'], ['A', 'I', 'C', 'E', 'H']]
goal = [['E', 'G', 'H'], ['A', 'D', 'I'], ['B', 'C', 'F']]

# Find and print solution
solution = find_solution(initial, goal)
print('<<<')
print('\n'.join(solution))
print('>>>')