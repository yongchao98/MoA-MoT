from heapq import heappop, heappush

def heuristic(state, goal_state):
    # Improved heuristic: count the number of blocks out of place
    score = 0
    for stack, goal_stack in zip(state, goal_state):
        for block, goal_block in zip(stack, goal_stack):
            if block != goal_block:
                score += 1
    return score

def blocks_world_optimized(initial_state, goal_state):
    # Priority queue for best-first search: (priority, current_state, moves)
    queue = [(heuristic(initial_state, goal_state), initial_state, [])]
    visited = set()
    visited.add(initial_state)

    while queue:
        _, current_state, moves = heappop(queue)

        # Check if the current state is the goal state
        if current_state == goal_state:
            return moves

        # Generate all possible valid moves
        for i, stack in enumerate(current_state):
            if stack:  # If the stack is not empty
                block = stack[-1]  # Top block
                for j in range(len(current_state)):
                    if i != j:  # Move to a different stack
                        new_state = list(map(list, current_state))
                        new_state[i].pop()
                        new_state[j].append(block)
                        new_state_tuple = tuple(map(tuple, new_state))

                        if new_state_tuple not in visited:
                            visited.add(new_state_tuple)
                            priority = heuristic(new_state_tuple, goal_state)
                            heappush(queue, (priority, new_state_tuple, moves + [f"Move {block} from {i+1} to {j+1}"]))

# Initial and goal states
initial_state = (('G', 'A', 'F', 'E', 'C'), ('J', 'H', 'B'), ('I', 'D'))
goal_state = (('A', 'B', 'D', 'F', 'I'), ('C', 'H'), ('E', 'G', 'J'))

# Find the sequence of moves
solution_moves = blocks_world_optimized(initial_state, goal_state)

# Print the solution
print("<<<" + "\n".join(solution_moves) + ">>>")