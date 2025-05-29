from heapq import heappop, heappush

def heuristic(state, goal_state):
    # Count the number of blocks not in their goal position
    misplaced = 0
    for stack, goal_stack in zip(state, goal_state):
        misplaced += sum(1 for block, goal_block in zip(stack, goal_stack) if block != goal_block)
    return misplaced

def blocksworld_astar(initial_state, goal_state):
    # Define a function to generate possible moves
    def generate_moves(state):
        moves = []
        for i, stack in enumerate(state):
            if stack:  # If the stack is not empty
                block = stack[-1]  # Get the top block
                for j in range(len(state)):
                    if i != j:  # Move to a different stack
                        new_state = list(map(list, state))
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