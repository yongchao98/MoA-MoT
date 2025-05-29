from collections import deque

def is_goal(state):
    # Correctly define the goal state
    return state == (['A', 'C', 'D'], ['B'], ['E'], [])

def get_valid_moves(state):
    stacks = state
    moves = []
    for i, stack in enumerate(stacks):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(stacks)):
                if i != j:  # Move to a different stack
                    new_stacks = [list(s) for s in stacks]  # Deep copy
                    new_stacks[j].append(new_stacks[i].pop())
                    move_description = f"Move {block} from stack{i+1} to stack{j+1}"
                    moves.append((move_description, tuple(map(tuple, new_stacks))))
    return moves

def bfs(initial_state):
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)

    while queue:
        current_state, path = queue.popleft()

        # Debugging: Print the current state and path
        print(f"Current state: {current_state}, Path: {path}")

        if is_goal(current_state):
            return path

        for move, new_state in get_valid_moves(current_state):
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + [move]))

    return None  # Return None if no solution is found

initial_state = (('D', 'C', 'E'), ('A',), (), ('B',))
solution = bfs(initial_state)

if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")