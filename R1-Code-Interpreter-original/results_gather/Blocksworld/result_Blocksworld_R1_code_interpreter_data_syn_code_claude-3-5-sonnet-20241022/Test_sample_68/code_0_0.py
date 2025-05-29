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
    new_stacks = copy.deepcopy(stacks)
    from_stack, to_stack, block = move
    new_stacks[from_stack].pop()
    new_stacks[to_stack].append(block)
    return new_stacks

def state_to_string(stacks):
    return ';'.join(','.join(stack) for stack in stacks)

def find_solution():
    initial = [['H', 'D'], ['E', 'J', 'A', 'C', 'F'], ['B', 'I', 'G', 'K']]
    goal = [['B', 'E', 'I', 'J', 'K'], ['A', 'C', 'D', 'F'], ['G', 'H']]
    
    queue = deque([(initial, [])])
    visited = {state_to_string(initial)}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if is_goal_state(current_state, goal):
            return moves
        
        for move in get_valid_moves(current_state):
            new_state = apply_move(current_state, move)
            state_str = state_to_string(new_state)
            
            if state_str not in visited:
                visited.add(state_str)
                from_stack, to_stack, block = move
                move_str = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                queue.append((new_state, moves + [move_str]))
    
    return None

solution = find_solution()
if solution:
    print("<<<")
    print('\n'.join(solution))
    print(">>>")
else:
    print("No solution found")