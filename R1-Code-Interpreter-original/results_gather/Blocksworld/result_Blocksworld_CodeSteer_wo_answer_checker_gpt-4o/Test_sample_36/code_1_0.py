from collections import deque

def bfs_with_heuristic(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = set()

    while queue:
        current_state, path = queue.popleft()
        state_tuple = tuple(map(tuple, current_state))

        if state_tuple in visited:
            continue

        visited.add(state_tuple)

        if current_state == goal_state:
            return path

        for i in range(len(current_state)):
            if not current_state[i]:
                continue
            block = current_state[i][-1]
            for j in range(len(current_state)):
                if i != j:
                    # Make a move
                    new_state = [stack[:] for stack in current_state]
                    new_state[i].pop()
                    new_state[j].append(block)
                    new_state_tuple = tuple(map(tuple, new_state))
                    if new_state_tuple not in visited:
                        move_description = f"Move {block} from stack{i+1} to stack{j+1}"
                        # Heuristic: prioritize moves that place blocks closer to their goal position
                        if block in goal_state[j] and new_state[j] == goal_state[j][:len(new_state[j])]:
                            queue.appendleft((new_state, path + [move_description]))
                        else:
                            queue.append((new_state, path + [move_description]))

    return None

initial_state = [['D', 'E'], ['G', 'I', 'F'], ['C', 'B', 'H', 'A']]
goal_state = [['B', 'C', 'E', 'F', 'I'], ['G', 'H'], ['A', 'D']]
solution = bfs_with_heuristic(initial_state, goal_state)

if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")