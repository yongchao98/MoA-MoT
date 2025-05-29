from collections import deque
import copy

def get_top_blocks(stacks):
    return [stack[-1] if stack else None for stack in stacks]

def is_goal_state(current, goal):
    return all(current[i] == goal[i] for i in range(len(current)))

def get_valid_moves(stacks):
    moves = []
    top_blocks = get_top_blocks(stacks)
    for i in range(len(stacks)):
        if not stacks[i]:
            continue
        for j in range(len(stacks)):
            if i != j:
                moves.append((i, j, stacks[i][-1]))
    return moves

def apply_move(stacks, move):
    from_stack, to_stack, block = move
    new_stacks = copy.deepcopy(stacks)
    new_stacks[from_stack].pop()
    new_stacks[to_stack].append(block)
    return new_stacks

def solve_blocksworld():
    initial = [['F', 'A', 'G'], ['E', 'B'], ['C', 'D']]
    goal = [['A', 'B', 'D'], ['G'], ['C', 'E', 'F']]
    
    queue = deque([(initial, [])])
    visited = set()
    
    while queue:
        current_stacks, path = queue.popleft()
        current_state_tuple = tuple(tuple(stack) for stack in current_stacks)
        
        if current_state_tuple in visited:
            continue
        
        visited.add(current_state_tuple)
        
        if is_goal_state(current_stacks, goal):
            return path
        
        for move in get_valid_moves(current_stacks):
            from_stack, to_stack, block = move
            new_stacks = apply_move(current_stacks, move)
            new_path = path + [f"Move {block} from {from_stack + 1} to {to_stack + 1}"]
            queue.append((new_stacks, new_path))
    
    return None

solution = solve_blocksworld()
if solution:
    print('\n'.join(solution))