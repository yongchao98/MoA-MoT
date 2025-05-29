def blocks_world_dfs(initial_state, goal_state):
    from collections import deque

    def is_goal(state):
        return state == goal_state

    def get_possible_moves(state):
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

    stack = deque([(initial_state, [])])
    visited = set()

    while stack:
        current_state, path = stack.pop()
        if is_goal(current_state):
            return path
        visited.add(current_state)
        for new_state, move in get_possible_moves(current_state):
            if new_state not in visited:
                stack.append((new_state, path + [move]))

initial_state = (('B', 'G', 'F'), ('A', 'E', 'C'), ('D',))
goal_state = (('D', 'F', 'G'), ('E',), ('A', 'B', 'C'))

solution = blocks_world_dfs(initial_state, goal_state)
print("<<<" + "\n".join(solution) + ">>>")