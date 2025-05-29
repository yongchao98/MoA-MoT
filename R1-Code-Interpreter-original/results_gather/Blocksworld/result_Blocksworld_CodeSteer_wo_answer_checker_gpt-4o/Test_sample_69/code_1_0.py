from heapq import heappop, heappush

def heuristic(state, goal_state):
    # Count the number of blocks not in their goal position or order
    misplaced = 0
    for stack, goal_stack in zip(state, goal_state):
        for block, goal_block in zip(stack, goal_stack):
            if block != goal_block:
                misplaced += 1
    return misplaced

def blocksworld_astar(initial_state, goal_state):
    # Priority queue for A* search
    open_set = []
    heappush(open_set, (0, initial_state, []))
    visited = set()
    visited.add(initial_state)

    while open_set:
        _, current_state, path = heappop(open_set)

        # Check if the current state is the goal state
        if current_state == goal_state:
            return path

        # Generate all possible moves
        for i, stack in enumerate(current_state):
            if stack:  # If the stack is not empty
                block = stack[-1]  # Get the top block
                for j in range(len(current_state)):
                    if i != j:  # Move to a different stack
                        new_state = list(map(list, current_state))
                        new_state[i].pop()
                        new_state[j].append(block)
                        new_state_tuple = tuple(map(tuple, new_state))

                        if new_state_tuple not in visited:
                            visited.add(new_state_tuple)
                            new_path = path + [f"Move {block} from {i+1} to {j+1}"]
                            cost = len(new_path) + heuristic(new_state_tuple, goal_state)
                            heappush(open_set, (cost, new_state_tuple, new_path))

# Define the initial and goal states
initial_state = (('H', 'I', 'B', 'E', 'A'), ('J',), ('G',), ('F', 'D', 'K', 'C'))
goal_state = (('H',), ('B', 'C', 'F', 'G', 'I'), ('A', 'D', 'E', 'J', 'K'))

# Find the sequence of moves
moves = blocksworld_astar(initial_state, goal_state)

# Print the moves
print("<<<" + "\n".join(moves) + ">>>")