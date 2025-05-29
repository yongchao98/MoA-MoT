from heapq import heappop, heappush

def heuristic(state, goal_state):
    # Simple heuristic: count the number of blocks not in their goal position
    misplaced = 0
    for stack, goal_stack in zip(state, goal_state):
        misplaced += sum(1 for block, goal_block in zip(stack, goal_stack) if block != goal_block)
    return misplaced

def a_star_blocks_world(initial_state, goal_state):
    # Priority queue for A*: each element is a tuple (cost, current_state, moves)
    open_set = []
    heappush(open_set, (0, initial_state, []))
    visited = set()
    visited.add(initial_state)

    while open_set:
        _, current_state, moves = heappop(open_set)

        # Check if the current state is the goal state
        if current_state == goal_state:
            return moves

        # Generate all possible moves
        for i, stack in enumerate(current_state):
            if stack:  # If the stack is not empty
                block = stack[-1]  # Get the top block
                for j in range(len(current_state)):
                    if i != j:  # Move to a different stack
                        new_state = list(map(list, current_state))
                        new_state[i].pop()  # Remove block from current stack
                        new_state[j].append(block)  # Add block to new stack
                        new_state_tuple = tuple(map(tuple, new_state))
                        if new_state_tuple not in visited:
                            visited.add(new_state_tuple)
                            new_cost = len(moves) + 1 + heuristic(new_state_tuple, goal_state)
                            heappush(open_set, (new_cost, new_state_tuple, moves + [f"Move {block} from {i+1} to {j+1}"]))

# Initial and goal states
initial_state = ((), ('E', 'B', 'I'), ('F', 'D'), ('C', 'H', 'G', 'A'))
goal_state = (('A', 'B', 'D', 'H'), ('I',), ('C', 'E', 'F', 'G'))

# Find the sequence of moves
moves = a_star_blocks_world(initial_state, goal_state)

# Print the moves in the required format
print("<<<" + "\n".join(moves) + ">>>")