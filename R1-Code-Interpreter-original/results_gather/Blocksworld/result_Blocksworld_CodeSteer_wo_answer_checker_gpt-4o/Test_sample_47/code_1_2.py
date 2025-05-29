def is_goal_state(stacks, goal):
    return stacks == goal

def move_block(stacks, from_stack, to_stack):
    # Make a deep copy of stacks to avoid modifying the original
    new_stacks = [stack[:] for stack in stacks]
    block = new_stacks[from_stack].pop()
    new_stacks[to_stack].append(block)
    return new_stacks

def dfs(stacks, goal, path):
    if is_goal_state(stacks, goal):
        return path

    for i in range(len(stacks)):
        if not stacks[i]:  # Skip empty stacks
            continue
        for j in range(len(stacks)):
            if i != j:  # Ensure we are not moving within the same stack
                new_stacks = move_block(stacks, i, j)
                new_path = path + [f"Move {stacks[i][-1]} from stack{i+1} to stack{j+1}"]
                result = dfs(new_stacks, goal, new_path)
                if result:
                    return result
    return None

# Initial state
initial_stacks = [
    ['D', 'G', 'E', 'C'],
    ['A'],
    ['F'],
    ['B']
]

# Goal state
goal_stacks = [
    ['D', 'G'],
    ['A', 'C'],
    ['B', 'E', 'F'],
    []
]

# Find the sequence of moves
solution = dfs(initial_stacks, goal_stacks, [])

# Print the solution
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")