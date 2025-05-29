from collections import deque
import copy

def get_top_blocks(state):
    tops = []
    for i, stack in enumerate(state):
        if stack:
            tops.append((stack[-1], i))
    return tops

def is_goal_state(current, goal):
    return current == goal

def get_valid_moves(state):
    moves = []
    top_blocks = get_top_blocks(state)
    
    # Try moving each top block to each other stack
    for block, source in top_blocks:
        for dest in range(len(state)):
            if dest != source:
                moves.append((block, source, dest))
    
    return moves

def apply_move(state, move):
    block, source, dest = move
    new_state = copy.deepcopy(state)
    block = new_state[source].pop()
    new_state[dest].append(block)
    return new_state

def state_to_string(state):
    return ';'.join(','.join(stack) for stack in state)

def find_solution(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = {state_to_string(initial_state)}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if is_goal_state(current_state, goal_state):
            return moves
        
        for move in get_valid_moves(current_state):
            new_state = apply_move(current_state, move)
            state_str = state_to_string(new_state)
            
            if state_str not in visited:
                visited.add(state_str)
                block, source, dest = move
                move_str = f"Move {block} from {source+1} to {dest+1}"
                queue.append((new_state, moves + [move_str]))
    
    return None

# Initial and goal states
initial = [
    ['D', 'E', 'H', 'G'],
    ['A', 'F', 'I'],
    ['B', 'C']
]

goal = [
    ['G', 'I'],
    ['A', 'C', 'D', 'F', 'H'],
    ['B', 'E']
]

# Find solution
solution = find_solution(initial, goal)

if solution:
    print("<<<")
    print('\n'.join(solution))
    print(">>>")
else:
    print("No solution found")