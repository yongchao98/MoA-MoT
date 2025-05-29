from heapq import heappop, heappush
from concurrent.futures import ThreadPoolExecutor
import itertools

def heuristic(state, goal):
    # Combine misplaced blocks and Manhattan distance
    misplaced = 0
    manhattan_distance = 0
    goal_positions = {block: (i, j) for i, stack in enumerate(goal) for j, block in enumerate(stack)}
    for i, stack in enumerate(state):
        for j, block in enumerate(stack):
            if block in goal_positions:
                goal_i, goal_j = goal_positions[block]
                if i != goal_i or j != goal_j:
                    misplaced += 1
                    manhattan_distance += abs(i - goal_i) + abs(j - goal_j)
    return misplaced + manhattan_distance

def get_successors(state):
    successors = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j in range(len(state)):
                if i != j:  # Move block to a different stack
                    new_state = list(map(list, state))
                    new_state[i].pop()
                    new_state[j].append(block)
                    successors.append((tuple(map(tuple, new_state)), f"Move {block} from {i+1} to {j+1}"))
    return successors

def a_star(initial, goal):
    open_set = []
    heappush(open_set, (0, initial, []))
    visited = set()
    best_cost = float('inf')

    while open_set:
        _, current, path = heappop(open_set)

        if current == goal:
            return path

        if current in visited:
            continue

        visited.add(current)

        for successor, move in get_successors(current):
            if successor not in visited:
                new_path = path + [move]
                cost = len(new_path) + heuristic(successor, goal)
                if cost < best_cost:
                    best_cost = cost
                    heappush(open_set, (cost, successor, new_path))

    return None

def parallel_a_star(initial, goal):
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(a_star, initial, goal)]
        for future in itertools.islice(futures, 1):
            solution = future.result()
            if solution:
                return solution
    return None

initial_state = (('F', 'A'), ('K', 'D', 'E', 'B', 'I', 'C'), ('H', 'J', 'G'))
goal_state = (('G',), ('B', 'D'), ('A', 'C', 'E', 'F', 'H', 'I', 'J', 'K'))

solution = parallel_a_star(initial_state, goal_state)

if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")