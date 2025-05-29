from collections import deque

def is_goal_state(state, goal_state):
    return state == goal_state

def get_possible_moves(state):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    new_state = [list(s) for s in state]  # Deep copy of the state
                    new_state[i].pop()  # Remove the block from the current stack
                    new_state[j].append(block)  # Add the block to the new stack
                    moves.append((new_state, f"Move {block} from {i+1} to {j+1}"))
    return moves

def bfs(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = set()

    while queue:
        current_state, path = queue.popleft()
        state_tuple = tuple(tuple(stack) for stack in current_state)

        if state_tuple in visited:
            continue

        visited.add(state_tuple)

        if is_goal_state(current_state, goal_state):
            return path

        for new_state, move in get_possible_moves(current_state):
            # Debugging output to trace the sequence of moves and states
            print(f"Current State: {current_state}")
            print(f"Move: {move}")
            print(f"New State: {new_state}\n")
            queue.append((new_state, path + [move]))

    return None

initial_state = [['A', 'C'], ['E', 'D', 'B', 'H', 'F'], ['G']]
goal_state = [['A', 'G', 'H'], ['F'], ['B', 'C', 'D', 'E']]

solution = bfs(initial_state, goal_state)

if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")