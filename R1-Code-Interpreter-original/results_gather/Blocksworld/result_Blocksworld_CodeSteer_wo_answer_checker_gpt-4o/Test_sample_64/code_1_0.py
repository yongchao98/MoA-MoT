from heapq import heappop, heappush

def heuristic(state, goal_state):
    # Calculate the sum of the Manhattan distances of each block from its goal position
    distance = 0
    goal_positions = {}
    for stack_index, stack in enumerate(goal_state):
        for block_index, block in enumerate(stack):
            goal_positions[block] = (stack_index, block_index)
    
    for stack_index, stack in enumerate(state):
        for block_index, block in enumerate(stack):
            if block in goal_positions:
                goal_stack_index, goal_block_index = goal_positions[block]
                distance += abs(stack_index - goal_stack_index) + abs(block_index - goal_block_index)
    
    return distance

def blocksworld_astar(initial_state, goal_state):
    # Define a function to generate possible moves
    def generate_moves(state):
        moves = []
        for i, stack in enumerate(state):
            if stack:  # If the stack is not empty
                block = stack[-1]  # Get the top block
                for j in range(len(state)):
                    if i != j:  # Move to a different stack
                        new_state = [list(s) for s in state]
                        new_state[i].pop()
                        new_state[j].append(block)
                        moves.append((tuple(map(tuple, new_state)), f"Move {block} from stack{i+1} to stack{j+1}"))
        return moves

    # Initialize the A* priority queue
    open_set = []
    heappush(open_set, (0, initial_state, []))
    visited = set()
    visited.add(initial_state)

    while open_set:
        _, current_state, path = heappop(open_set)

        # Check if the current state is the goal state
        if current_state == goal_state:
            return path

        # Generate and explore new states
        for new_state, move in generate_moves(current_state):
            if new_state not in visited:
                visited.add(new_state)
                cost = len(path) + 1 + heuristic(new_state, goal_state)
                heappush(open_set, (cost, new_state, path + [move]))

# Define the initial and goal states
initial_state = (('G', 'C'), ('D', 'A', 'E', 'B'), ('I', 'H', 'F'))
goal_state = (('B', 'I'), ('E',), ('A', 'C', 'D', 'F', 'G', 'H'))

# Find the sequence of moves
moves = blocksworld_astar(initial_state, goal_state)

# Print the moves
print("<<<" + "\n".join(moves) + ">>>")