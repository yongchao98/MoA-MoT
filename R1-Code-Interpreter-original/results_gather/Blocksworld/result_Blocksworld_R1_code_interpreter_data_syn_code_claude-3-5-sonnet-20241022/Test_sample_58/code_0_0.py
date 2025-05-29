from collections import deque
import copy

def get_top_block(stack):
    return stack[-1] if stack else None

def is_goal_state(current, goal):
    return all(current[i] == goal[i] for i in range(len(current)))

def get_valid_moves(stacks):
    moves = []
    for i, source in enumerate(stacks):
        if source:  # if stack is not empty
            top_block = source[-1]
            for j, dest in enumerate(stacks):
                if i != j:  # can't move to same stack
                    moves.append((i, j, top_block))
    return moves

def apply_move(stacks, move):
    source, dest, block = move
    new_stacks = copy.deepcopy(stacks)
    new_stacks[source].pop()
    new_stacks[dest].append(block)
    return new_stacks

def state_to_string(stacks):
    return ';'.join(','.join(stack) for stack in stacks)

def find_solution():
    initial = [['C', 'D', 'B'], ['E'], ['F', 'A']]
    goal = [['A', 'C', 'F'], ['B'], ['D', 'E']]
    
    queue = deque([(initial, [])])
    visited = {state_to_string(initial)}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if is_goal_state(current_state, goal):
            return moves
        
        for move in get_valid_moves(current_state):
            source, dest, block = move
            new_state = apply_move(current_state, move)
            state_str = state_to_string(new_state)
            
            if state_str not in visited:
                visited.add(state_str)
                move_str = f"Move {block} from {source+1} to {dest+1}"
                queue.append((new_state, moves + [move_str]))
    
    return None

solution = find_solution()
if solution:
    print('\n'.join(solution))