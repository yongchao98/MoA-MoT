from collections import deque

def blocksworld_bfs(initial_state, goal_state):
    # Helper function to generate new states
    def generate_new_states(state):
        new_states = []
        for i, stack in enumerate(state):
            if stack:  # If the stack is not empty
                block = stack[-1]  # Get the top block
                for j in range(len(state)):
                    if i != j:  # Move to a different stack
                        new_state = list(map(list, state))
                        new_state[i].pop()  # Remove block from current stack
                        new_state[j].append(block)  # Add block to new stack
                        new_states.append((tuple(map(tuple, new_state)), f"Move {block} from {i+1} to {j+1}"))
        return new_states

    # Heuristic function to estimate distance to goal
    def heuristic(state):
        distance = 0
        for i, stack in enumerate(state):
            for j, block in enumerate(stack):
                if j < len(goal_state[i]) and block != goal_state[i][j]:
                    distance += 1
        return distance

    # BFS setup with priority queue
    queue = deque([(initial_state, [], 0)])
    visited = set()
    visited.add(initial_state)

    while queue:
        # Sort queue based on heuristic value
        queue = deque(sorted(queue, key=lambda x: x[2]))

        current_state, moves, _ = queue.popleft()

        if current_state == goal_state:
            return moves

        for new_state, move in generate_new_states(current_state):
            if new_state not in visited:
                visited.add(new_state)
                new_moves = moves + [move]
                queue.append((new_state, new_moves, heuristic(new_state)))

    return None

# Initial and goal states
initial_state = (('G', 'I', 'D'), ('E', 'A', 'K'), ('C', 'B', 'H', 'J', 'F'))
goal_state = (('A', 'D', 'F'), ('B', 'H', 'I', 'J', 'K'), ('C', 'E', 'G'))

# Find the sequence of moves
moves = blocksworld_bfs(initial_state, goal_state)

# Output the moves
if moves:
    print("<<<" + "\n".join(moves) + ">>>")
else:
    print("No solution found.")