import heapq

def heuristic(state, goal_state):
    """Estimate the cost to reach the goal state."""
    cost = 0
    for stack, goal_stack in zip(state, goal_state):
        for block, goal_block in zip(stack, goal_stack):
            if block != goal_block:
                cost += 1
    return cost

def move_block(state, from_stack, to_stack):
    """Move the top block from one stack to another."""
    new_state = [list(stack) for stack in state]
    block = new_state[from_stack].pop()
    new_state[to_stack].append(block)
    return tuple(tuple(stack) for stack in new_state), block

def is_goal_state(state, goal_state):
    """Check if the current state matches the goal state."""
    return state == goal_state

def solve_blocks(initial_state, goal_state):
    """Solve the blocksworld problem using A* search."""
    initial_state = tuple(tuple(stack) for stack in initial_state)
    goal_state = tuple(tuple(stack) for stack in goal_state)
    
    # Priority queue for A* search
    open_set = []
    heapq.heappush(open_set, (0, initial_state, []))
    
    # Set to track visited states
    visited = set()
    
    while open_set:
        _, current_state, moves = heapq.heappop(open_set)
        
        if is_goal_state(current_state, goal_state):
            return moves
        
        if current_state in visited:
            continue
        
        visited.add(current_state)
        
        for from_stack in range(len(current_state)):
            if not current_state[from_stack]:
                continue
            for to_stack in range(len(current_state)):
                if from_stack == to_stack:
                    continue
                # Generate new state
                new_state, block = move_block(current_state, from_stack, to_stack)
                new_moves = moves + [f"Move {block} from {from_stack + 1} to {to_stack + 1}"]
                
                if new_state not in visited:
                    cost = len(new_moves) + heuristic(new_state, goal_state)
                    heapq.heappush(open_set, (cost, new_state, new_moves))
    
    return None

# Initial and goal states
initial_state = [['G', 'E', 'B'], ['A'], ['C', 'D', 'H', 'F']]
goal_state = [['D', 'F', 'G'], ['A', 'C', 'E', 'H'], ['B']]

# Solve the problem
moves = solve_blocks(initial_state, goal_state)
if moves:
    print("<<<")
    for move in moves:
        print(move)
    print(">>>")
else:
    print("No solution found.")