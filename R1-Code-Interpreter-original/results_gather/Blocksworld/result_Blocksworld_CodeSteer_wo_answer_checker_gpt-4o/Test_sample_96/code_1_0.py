import heapq
from multiprocessing import Pool, cpu_count

def heuristic(state, goal_state):
    # Calculate a heuristic based on the number of misplaced blocks and their distance from the goal
    misplaced = 0
    for stack_index, (stack, goal_stack) in enumerate(zip(state, goal_state)):
        for block_index, block in enumerate(stack):
            if block_index < len(goal_stack) and block != goal_stack[block_index]:
                misplaced += 1
            elif block_index >= len(goal_stack):
                misplaced += 1
    return misplaced

def generate_moves(current_state):
    moves = []
    for i, stack in enumerate(current_state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(current_state)):
                if i != j:  # Move to a different stack
                    new_state = list(map(list, current_state))
                    new_state[i].pop()  # Remove block from current stack
                    new_state[j].append(block)  # Add block to new stack
                    new_state_tuple = tuple(map(tuple, new_state))
                    move = f"Move {block} from stack{i+1} to stack{j+1}"
                    moves.append((new_state_tuple, move))
    return moves

def blocks_world_astar(initial_state, goal_state):
    # Priority queue for A*: stores (estimated_cost, current_cost, current_state, path_to_state)
    priority_queue = []
    heapq.heappush(priority_queue, (0, 0, initial_state, []))
    visited = set()

    while priority_queue:
        estimated_cost, current_cost, current_state, path = heapq.heappop(priority_queue)

        # Check if the current state is the goal state
        if current_state == goal_state:
            return path

        # Mark the current state as visited
        visited.add(current_state)

        # Generate all possible moves
        for new_state_tuple, move in generate_moves(current_state):
            if new_state_tuple not in visited:
                new_cost = current_cost + 1
                estimated_total_cost = new_cost + heuristic(new_state_tuple, goal_state)
                heapq.heappush(priority_queue, (estimated_total_cost, new_cost, new_state_tuple, path + [move]))

    return None  # If no solution is found

# Initial and goal states
initial_state = (('E', 'B', 'J', 'I'), ('A', 'F', 'K', 'C'), ('D', 'G', 'H'))
goal_state = (('A', 'C', 'F', 'J', 'K'), ('D', 'G'), ('B', 'E', 'H', 'I'))

# Find the sequence of moves
solution = blocks_world_astar(initial_state, goal_state)

# Print the solution
if solution:
    print("<<<")
    for move in solution:
        print(move)
    print(">>>")
else:
    print("No solution found.")