import heapq

def heuristic(state, goal):
    # Count the number of blocks out of place
    misplaced = 0
    for stack_id, stack in enumerate(state):
        goal_stack = goal[stack_id]
        for i, block in enumerate(stack):
            if i >= len(goal_stack) or block != goal_stack[i]:
                misplaced += 1
    return misplaced

def get_possible_moves(state):
    moves = []
    for src in range(len(state)):
        if state[src]:  # If the stack is not empty
            block = state[src][-1]
            for dest in range(len(state)):
                if src != dest:
                    moves.append((block, src, dest))
    return moves

def apply_move(state, move):
    block, src, dest = move
    new_state = [list(stack) for stack in state]
    new_state[src].pop()
    new_state[dest].append(block)
    return tuple(tuple(stack) for stack in new_state)

def a_star_search(initial_state, goal_state):
    open_set = []
    heapq.heappush(open_set, (0, initial_state, []))
    visited = set()

    while open_set:
        _, current_state, path = heapq.heappop(open_set)

        if current_state == goal_state:
            return path

        if current_state in visited:
            continue
        visited.add(current_state)

        for move in get_possible_moves(current_state):
            new_state = apply_move(current_state, move)
            new_path = path + [move]
            cost = len(new_path) + heuristic(new_state, goal_state)
            heapq.heappush(open_set, (cost, new_state, new_path))

    return None

def blocksworld_solution():
    initial_state = (
        ('G',),
        ('C', 'B'),
        ('E',),
        ('A', 'F', 'D')
    )

    goal_state = (
        ('A', 'C'),
        ('B', 'D', 'E', 'G'),
        ('F',),
        ()
    )

    solution = a_star_search(initial_state, goal_state)
    if solution:
        for block, src, dest in solution:
            print(f"Move {block} from {src + 1} to {dest + 1}")
    else:
        print("No solution found")

blocksworld_solution()