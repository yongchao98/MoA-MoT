import heapq
import time

def is_goal_state(stacks, goal):
    return stacks == goal

def get_valid_moves(stacks):
    moves = []
    for i, stack in enumerate(stacks):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j in range(len(stacks)):
                if i != j:  # Move to a different stack
                    moves.append((block, i, j))
    return moves

def apply_move(stacks, move):
    block, src, dest = move
    new_stacks = [stack[:] for stack in stacks]
    new_stacks[src].pop()
    new_stacks[dest].append(block)
    return new_stacks

def heuristic(stacks, goal):
    # Improved heuristic: sum of distances each block is from its goal position
    distance = 0
    goal_positions = {block: (i, j) for i, stack in enumerate(goal) for j, block in enumerate(stack)}
    for i, stack in enumerate(stacks):
        for j, block in enumerate(stack):
            if block in goal_positions:
                goal_i, goal_j = goal_positions[block]
                distance += abs(i - goal_i) + abs(j - goal_j)
    return distance

def a_star(initial_state, goal_state, max_iterations=20000):
    start_time = time.time()
    open_set = []
    heapq.heappush(open_set, (0, initial_state, []))
    visited = set()
    visited.add(tuple(tuple(stack) for stack in initial_state))

    iterations = 0
    while open_set and iterations < max_iterations:
        _, current_stacks, path = heapq.heappop(open_set)

        if is_goal_state(current_stacks, goal_state):
            return path

        for move in get_valid_moves(current_stacks):
            new_stacks = apply_move(current_stacks, move)
            new_state = tuple(tuple(stack) for stack in new_stacks)

            if new_state not in visited:
                visited.add(new_state)
                cost = len(path) + 1 + heuristic(new_stacks, goal_state)
                heapq.heappush(open_set, (cost, new_stacks, path + [move]))

        iterations += 1
        if iterations % 1000 == 0:
            print(f"Iteration: {iterations}, Open set size: {len(open_set)}")

    return None

initial_state = [['G', 'A', 'H'], ['E'], ['I', 'C', 'B'], ['F', 'D']]
goal_state = [['A', 'B', 'C', 'D', 'G', 'H', 'I'], ['E'], ['F']]

solution = a_star(initial_state, goal_state)

if solution:
    print("<<<")
    for block, src, dest in solution:
        print(f"Move {block} from {src + 1} to {dest + 1}")
    print(">>>")
else:
    print("No solution found or timeout reached.")