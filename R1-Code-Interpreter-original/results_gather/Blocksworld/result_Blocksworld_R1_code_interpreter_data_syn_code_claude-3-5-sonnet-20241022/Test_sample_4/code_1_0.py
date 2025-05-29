from collections import deque
import copy

def get_top_blocks(stacks):
    return [stack[-1] if stack else None for stack in stacks]

def is_goal_state(current, goal):
    return all(current[i] == goal[i] for i in range(len(current)))

def get_valid_moves(stacks):
    moves = []
    for i in range(len(stacks)):
        if not stacks[i]:  # Skip empty stacks
            continue
        moving_block = stacks[i][-1]  # Can only move top block
        for j in range(len(stacks)):
            if i != j:
                moves.append((i, j, moving_block))
    return moves

def state_to_string(stacks):
    return ';'.join(','.join(stack) for stack in stacks)

def solve_blocksworld():
    initial = [['F', 'A', 'G'], ['E', 'B'], ['C', 'D']]
    goal = [['A', 'B', 'D'], ['G'], ['C', 'E', 'F']]
    
    queue = deque([(initial, [])])
    visited = set()
    
    while queue:
        current_stacks, path = queue.popleft()
        current_state = state_to_string(current_stacks)
        
        if current_state in visited:
            continue
            
        if is_goal_state(current_stacks, goal):
            return path
            
        visited.add(current_state)
        
        for move in get_valid_moves(current_stacks):
            from_stack, to_stack, block = move
            new_stacks = copy.deepcopy(current_stacks)
            # Only move if it's actually the top block
            if new_stacks[from_stack][-1] == block:
                block = new_stacks[from_stack].pop()
                new_stacks[to_stack].append(block)
                new_path = path + [f"Move {block} from {from_stack + 1} to {to_stack + 1}"]
                queue.append((new_stacks, new_path))

solution = solve_blocksworld()
if solution:
    print('\n'.join(solution))