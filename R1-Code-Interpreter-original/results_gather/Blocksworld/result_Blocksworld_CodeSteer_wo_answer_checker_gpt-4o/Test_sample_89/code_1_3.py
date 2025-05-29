import heapq

def heuristic(state, goal):
    # Simple heuristic: count the number of blocks out of place
    misplaced = 0
    for stack_id, stack in state.items():
        goal_stack = goal.get(stack_id, [])
        for i, block in enumerate(stack):
            if i >= len(goal_stack) or block != goal_stack[i]:
                misplaced += 1
    return misplaced

def get_possible_moves(state):
    moves = []
    for src in state:
        if state[src]:  # If the stack is not empty
            block = state[src][-1]
            for dest in state:
                if src != dest:
                    moves.append((block, src, dest))
    return moves

def apply_move(state, move):
    block, src, dest = move
    new_state = {k: v[:] for k, v in state.items()}
    new_state[src].pop()
    new_state[dest].append(block)
    return new_state

def a_star_search(initial_state, goal_state):
    open_set = []
    heapq.heappush(open_set, (0, initial_state, []))
    visited = set()

    while open_set:
        _, current_state, path = heapq.heappop(open_set)

        if current_state == goal_state:
            return path

        state_tuple = tuple(tuple(stack) for stack in current_state.values())
        if state_tuple in visited:
            continue
        visited.add(state_tuple)

        for move in get_possible_moves(current_state):
            new_state = apply_move(current_state, move)
            new_path = path + [move]
            cost = len(new_path) + heuristic(new_state, goal_state)
            heapq.heappush(open_set, (cost, new_state, new_path))

    return None

def blocksworld_solution():
    initial_state = {
        1: ['G'],
        2: ['C', 'B'],
        3: ['E'],
        4: ['A', 'F', 'D']
    }

    goal_state = {
        1: ['A', 'C'],
        2: ['B', 'D', 'E', 'G'],
        3: ['F'],
        4: []
    }

    solution = a_star_search(initial_state, goal_state)
    if solution:
        for block, src, dest in solution:
            print(f"Move {block} from {src} to {dest}")
    else:
        print("No solution found")

blocksworld_solution()