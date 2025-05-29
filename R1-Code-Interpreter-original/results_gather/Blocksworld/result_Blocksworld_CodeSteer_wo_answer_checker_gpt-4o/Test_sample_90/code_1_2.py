from collections import deque

def blocksworld_bfs(initial_state, goal_state):
    def is_goal(state):
        return state == goal_state

    def get_possible_moves(state):
        moves = []
        for i, stack in enumerate(state):
            if stack:  # If the stack is not empty
                block = stack[-1]  # Get the top block
                for j in range(len(state)):
                    if i != j:  # Move to a different stack
                        new_state = [list(s) for s in state]
                        new_state[i].pop()
                        new_state[j].append(block)
                        moves.append((tuple(map(tuple, new_state)), f"Move {block} from {i+1} to {j+1}"))
        return moves

    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)

    while queue:
        current_state, path = queue.popleft()

        if is_goal(current_state):
            return path

        for new_state, move in get_possible_moves(current_state):
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + [move]))

    return None

initial_state = (('D', 'E', 'G'), ('C', 'B', 'A'), ('F', 'H'))
goal_state = (('B', 'E', 'G'), ('C', 'D', 'H'), ('A', 'F'))

solution = blocksworld_bfs(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")