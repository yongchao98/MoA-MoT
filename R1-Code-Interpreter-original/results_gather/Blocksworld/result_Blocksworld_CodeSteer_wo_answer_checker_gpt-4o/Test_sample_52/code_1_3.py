import heapq

def heuristic(state, goal):
    # Calculate the heuristic based on misplaced blocks and their distances
    misplaced = 0
    for stack_index, (stack, goal_stack) in enumerate(zip(state, goal)):
        for block_index, block in enumerate(stack):
            if block_index < len(goal_stack) and block != goal_stack[block_index]:
                misplaced += 1
            elif block_index >= len(goal_stack):
                misplaced += 1
    return misplaced

def get_possible_moves(state):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    new_state = [list(s) for s in state]
                    new_state[i].pop()
                    new_state[j].append(block)
                    moves.append((tuple(tuple(s) for s in new_state), f"Move {block} from {i + 1} to {j + 1}"))
    return moves

def a_star(initial, goal):
    open_set = []
    heapq.heappush(open_set, (0, initial, []))
    closed_set = set()

    while open_set:
        _, current_state, path = heapq.heappop(open_set)

        if current_state == goal:
            return path

        if current_state in closed_set:
            continue

        closed_set.add(current_state)

        for next_state, move in get_possible_moves(current_state):
            if next_state not in closed_set:
                new_path = path + [move]
                cost = len(new_path) + heuristic(next_state, goal)
                heapq.heappush(open_set, (cost, next_state, new_path))

    return None

initial_state = (
    ('I', 'J', 'F'),
    ('B', 'G', 'E', 'C'),
    ('H', 'D', 'A')
)

goal_state = (
    ('B', 'C', 'H', 'I', 'J'),
    ('D',),
    ('A', 'E', 'F', 'G')
)

moves = a_star(initial_state, goal_state)
if moves:
    print("<<<" + "\n".join(moves) + ">>>")
else:
    print("No solution found.")