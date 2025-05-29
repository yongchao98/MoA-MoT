import heapq
import time

def is_goal_state(stacks, goal):
    return stacks == goal

def get_possible_moves(stacks):
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
    new_stacks = [list(stack) for stack in stacks]  # Deep copy
    new_stacks[src].pop()
    new_stacks[dest].append(block)
    return new_stacks

def heuristic(stacks, goal):
    # Simple heuristic: count the number of blocks out of place
    misplaced = 0
    for stack, goal_stack in zip(stacks, goal):
        misplaced += sum(1 for block, goal_block in zip(stack, goal_stack) if block != goal_block)
    return misplaced

def a_star_search(initial_state, goal_state, time_limit=10):
    start_time = time.time()
    priority_queue = []
    heapq.heappush(priority_queue, (0, initial_state, []))
    visited = set()

    while priority_queue:
        if time.time() - start_time > time_limit:
            print("Timeout reached, no solution found within time limit.")
            return None

        _, current_state, path = heapq.heappop(priority_queue)
        state_tuple = tuple(tuple(stack) for stack in current_state)

        if state_tuple in visited:
            continue

        visited.add(state_tuple)

        if is_goal_state(current_state, goal_state):
            return path

        for move in get_possible_moves(current_state):
            new_state = apply_move(current_state, move)
            new_state_tuple = tuple(tuple(stack) for stack in new_state)
            if new_state_tuple not in visited:
                cost = len(path) + 1 + heuristic(new_state, goal_state)
                heapq.heappush(priority_queue, (cost, new_state, path + [move]))

    return None

initial_state = [['B', 'A'], ['C', 'F', 'G'], ['D', 'I', 'E', 'J', 'H']]
goal_state = [['A', 'E', 'I'], ['B', 'C', 'J'], ['D', 'F', 'G', 'H']]

solution = a_star_search(initial_state, goal_state)

if solution:
    print("<<<")
    for block, src, dest in solution:
        print(f"Move {block} from stack{src+1} to stack{dest+1}")
    print(">>>")
else:
    print("No solution found.")