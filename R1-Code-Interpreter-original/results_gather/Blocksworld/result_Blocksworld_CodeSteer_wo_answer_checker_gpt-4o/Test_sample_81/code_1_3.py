from collections import deque
from multiprocessing import Pool, cpu_count

def is_goal_state(state, goal_state):
    return state == goal_state

def get_valid_moves(state):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    moves.append((block, i, j))
    return moves

def apply_move(state, move):
    block, src, dest = move
    new_state = [stack[:] for stack in state]  # Deep copy of the state
    new_state[src].pop()
    new_state[dest].append(block)
    return new_state

def heuristic(state, goal_state):
    # Simple heuristic: count the number of blocks in the correct stack
    score = 0
    for stack, goal_stack in zip(state, goal_state):
        score += sum(1 for block in stack if block in goal_stack)
    return score

def bfs(initial_state, goal_state, max_depth=20):
    queue = deque([(initial_state, [], 0)])
    visited = set()

    while queue:
        current_state, path, depth = queue.popleft()

        if is_goal_state(current_state, goal_state):
            return path

        if depth > max_depth:
            continue

        state_tuple = tuple(tuple(stack) for stack in current_state)
        if state_tuple in visited:
            continue

        visited.add(state_tuple)

        valid_moves = get_valid_moves(current_state)
        valid_moves.sort(key=lambda move: heuristic(apply_move(current_state, move), goal_state), reverse=True)

        for move in valid_moves:
            new_state = apply_move(current_state, move)
            queue.append((new_state, path + [move], depth + 1))

    return None

def parallel_bfs(initial_state, goal_state):
    with Pool(cpu_count()) as pool:
        results = pool.starmap(bfs, [(initial_state, goal_state) for _ in range(cpu_count())])
        for result in results:
            if result:
                return result
    return None

initial_state = [['D'], ['A', 'J'], ['H', 'F', 'C', 'E', 'I'], ['B', 'G']]
goal_state = [['F', 'I'], ['A', 'B', 'C', 'D', 'G'], ['E', 'H', 'J']]

solution = parallel_bfs(initial_state, goal_state)

if solution:
    print("<<<" + "\n".join(f"Move {block} from stack{src+1} to stack{dest+1}" for block, src, dest in solution) + ">>>")
else:
    print("No solution found.")