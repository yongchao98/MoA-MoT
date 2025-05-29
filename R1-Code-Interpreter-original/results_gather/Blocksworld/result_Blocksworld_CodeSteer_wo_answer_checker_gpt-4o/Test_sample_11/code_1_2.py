# Define the initial and goal states
initial_state = {
    'stack1': [],
    'stack2': ['E'],
    'stack3': ['B', 'A', 'F', 'J'],
    'stack4': ['C', 'G', 'D', 'H', 'I']
}

goal_state = {
    'stack1': ['C', 'D'],
    'stack2': ['E', 'F', 'I'],
    'stack3': ['A', 'B', 'G', 'H', 'J']
}

# Function to perform a move
def move_block(state, source, destination):
    block = state[source].pop()
    state[destination].append(block)
    return f"Move {block} from {source} to {destination}"

# Function to check if the current state matches the goal state
def is_goal_state(state, goal_state):
    return state == goal_state

# Function to perform DFS with backtracking
def dfs(state, goal_state, path, visited):
    if is_goal_state(state, goal_state):
        return path

    for source in state:
        if state[source]:  # If the stack is not empty
            block = state[source][-1]
            for destination in state:
                if source != destination:
                    # Perform the move
                    move = move_block(state, source, destination)
                    path.append(move)

                    # Check if the new state has been visited
                    state_tuple = tuple(tuple(state[stack]) for stack in state)
                    if state_tuple not in visited:
                        visited.add(state_tuple)
                        result = dfs(state, goal_state, path, visited)
                        if result:
                            return result

                    # Backtrack
                    path.pop()
                    state[destination].pop()
                    state[source].append(block)

    return None

# Initialize the path and visited set
path = []
visited = set()
state_tuple = tuple(tuple(initial_state[stack]) for stack in initial_state)
visited.add(state_tuple)

# Perform DFS
solution = dfs(initial_state, goal_state, path, visited)

# Print the solution
if solution:
    print("Solution found:")
    print("<<<")
    for move in solution:
        print(move)
    print(">>>")
else:
    print("No solution found.")