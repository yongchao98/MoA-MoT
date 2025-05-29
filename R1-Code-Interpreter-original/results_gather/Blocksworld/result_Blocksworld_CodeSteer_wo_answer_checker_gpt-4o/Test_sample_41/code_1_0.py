import heapq

def is_goal_state(stacks, goal_stacks):
    return stacks == goal_stacks

def heuristic(stacks, goal_stacks):
    # Count the number of blocks not in their goal position
    misplaced_blocks = 0
    for stack, goal_stack in zip(stacks, goal_stacks):
        for block, goal_block in zip(stack, goal_stack):
            if block != goal_block:
                misplaced_blocks += 1
    return misplaced_blocks

def get_possible_moves(stacks):
    moves = []
    for i, stack in enumerate(stacks):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j, target_stack in enumerate(stacks):
                if i != j:  # Don't move to the same stack
                    moves.append((block, i, j))
    return moves

def apply_move(stacks, move):
    block, source, destination = move
    new_stacks = [stack[:] for stack in stacks]
    new_stacks[source].pop()
    new_stacks[destination].append(block)
    return new_stacks

def a_star_search(initial_stacks, goal_stacks):
    open_set = []
    heapq.heappush(open_set, (0, initial_stacks, []))
    visited = set()
    initial_tuple = tuple(tuple(stack) for stack in initial_stacks)
    visited.add(initial_tuple)

    while open_set:
        _, current_stacks, path = heapq.heappop(open_set)

        if is_goal_state(current_stacks, goal_stacks):
            return path

        for move in get_possible_moves(current_stacks):
            new_stacks = apply_move(current_stacks, move)
            stacks_tuple = tuple(tuple(stack) for stack in new_stacks)
            if stacks_tuple not in visited:
                visited.add(stacks_tuple)
                new_path = path + [move]
                cost = len(new_path) + heuristic(new_stacks, goal_stacks)
                heapq.heappush(open_set, (cost, new_stacks, new_path))

    return None

initial_stacks = [
    ['I'],
    ['G', 'C', 'H', 'A', 'E'],
    ['F', 'B'],
    ['D', 'K', 'J']
]

goal_stacks = [
    ['A', 'D', 'F', 'I', 'J'],
    ['G', 'H'],
    ['B', 'C', 'E', 'K']
]

solution = a_star_search(initial_stacks, goal_stacks)

if solution:
    print("<<<")
    for move in solution:
        block, source, destination = move
        print(f"Move {block} from {source + 1} to {destination + 1}")
    print(">>>")
else:
    print("No solution found.")