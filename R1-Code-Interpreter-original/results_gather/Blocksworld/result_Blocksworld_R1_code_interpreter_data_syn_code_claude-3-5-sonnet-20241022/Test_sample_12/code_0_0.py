from collections import deque
import copy

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty and destination stack exists
    if not stacks[from_stack] or to_stack not in [0, 1, 2]:
        return False
    return True

def get_state_string(stacks):
    return "|".join(",".join(stack) for stack in stacks)

def get_moves(current, goal):
    # Initialize stacks
    initial_stacks = [[] for _ in range(3)]
    for i, stack in enumerate(current):
        initial_stacks[i].extend(list(stack))
    
    goal_stacks = [[] for _ in range(3)]
    for i, stack in enumerate(goal):
        goal_stacks[i].extend(list(stack))
    
    # BFS search
    queue = deque([(initial_stacks, [])])
    visited = {get_state_string(initial_stacks)}
    
    while queue:
        current_stacks, moves = queue.popleft()
        
        # Check if goal reached
        if current_stacks == goal_stacks:
            return moves
        
        # Try all possible moves
        for from_stack in range(3):
            for to_stack in range(3):
                if from_stack != to_stack and is_valid_move(current_stacks, from_stack, to_stack):
                    # Create new state
                    new_stacks = copy.deepcopy(current_stacks)
                    block = new_stacks[from_stack].pop()
                    new_stacks[to_stack].append(block)
                    
                    state_string = get_state_string(new_stacks)
                    if state_string not in visited:
                        visited.add(state_string)
                        new_moves = moves + [f"Move {block} from {from_stack + 1} to {to_stack + 1}"]
                        queue.append((new_stacks, new_moves))
    
    return None

# Initial and goal states
initial = [
    [],
    ['K', 'D', 'I', 'G', 'B', 'J'],
    ['H', 'C', 'E', 'A', 'F']
]

goal = [
    ['A', 'E', 'F'],
    ['B', 'G', 'K'],
    ['C', 'D', 'H', 'I', 'J']
]

# Find solution
solution = get_moves(initial, goal)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")