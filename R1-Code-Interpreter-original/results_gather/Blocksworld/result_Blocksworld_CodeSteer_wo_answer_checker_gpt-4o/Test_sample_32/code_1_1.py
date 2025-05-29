from collections import deque

def bfs_blocks_world(initial_state, goal_state):
    def is_goal(state):
        return state == goal_state

    def get_valid_moves(state):
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

    queue = deque([(initial_state, [])])
    visited = set()

    while queue:
        current_state, path = queue.popleft()
        if is_goal(current_state):
            return path

        if current_state in visited:
            continue

        visited.add(current_state)

        for new_state, move in get_valid_moves(current_state):
            if new_state not in visited:
                queue.append((new_state, path + [move]))

    return None

initial_state = (('G', 'D', 'B', 'C'), ('F', 'A', 'E'), ())
goal_state = (('A', 'B', 'D', 'E', 'G'), ('C',), ('F',))

solution = bfs_blocks_world(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")