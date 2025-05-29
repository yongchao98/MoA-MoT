from collections import deque

def is_goal_state(stacks, goal):
    return stacks == goal

def get_possible_moves(stacks):
    moves = []
    for i, stack in enumerate(stacks):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j in range(len(stacks)):
                if i != j:  # Move to a different stack
                    new_stacks = [list(s) for s in stacks]  # Deep copy
                    new_stacks[i].pop()
                    new_stacks[j].append(block)
                    moves.append((new_stacks, f"Move {block} from {i+1} to {j+1}"))
    return moves

def bfs(initial, goal):
    queue = deque([(initial, [])])
    visited = set()

    while queue:
        current_stacks, path = queue.popleft()
        if is_goal_state(current_stacks, goal):
            return path

        state_tuple = tuple(tuple(stack) for stack in current_stacks)
        if state_tuple in visited:
            continue
        visited.add(state_tuple)

        for new_stacks, move in get_possible_moves(current_stacks):
            queue.append((new_stacks, path + [move]))

    return None

initial_state = [['A', 'D'], [], ['F', 'B'], ['E', 'C']]
goal_state = [['E'], ['A', 'B', 'C', 'F'], ['D'], []]

solution = bfs(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")