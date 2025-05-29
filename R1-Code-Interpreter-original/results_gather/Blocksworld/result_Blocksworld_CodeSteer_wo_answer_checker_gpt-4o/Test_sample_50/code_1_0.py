import heapq

def is_goal_state(stacks, goal_state):
    return stacks == goal_state

def move_block(stacks, from_stack, to_stack):
    new_stacks = [list(stack) for stack in stacks]
    block = new_stacks[from_stack].pop()
    new_stacks[to_stack].append(block)
    return new_stacks

def stacks_to_tuple(stacks):
    return tuple(tuple(stack) for stack in stacks)

def heuristic(stacks, goal_state):
    # Simple heuristic: count the number of misplaced blocks
    misplaced = 0
    for stack, goal_stack in zip(stacks, goal_state):
        misplaced += sum(1 for block, goal_block in zip(stack, goal_stack) if block != goal_block)
    return misplaced

def a_star_search(initial_state, goal_state):
    open_set = []
    heapq.heappush(open_set, (0, initial_state, []))
    visited = set()

    while open_set:
        _, current_stacks, path = heapq.heappop(open_set)
        current_tuple = stacks_to_tuple(current_stacks)

        if current_tuple in visited:
            continue
        visited.add(current_tuple)

        if is_goal_state(current_stacks, goal_state):
            return path

        for from_stack in range(len(current_stacks)):
            if not current_stacks[from_stack]:
                continue
            for to_stack in range(len(current_stacks)):
                if from_stack == to_stack:
                    continue
                new_stacks = move_block(current_stacks, from_stack, to_stack)
                new_path = path + [f"Move {current_stacks[from_stack][-1]} from stack{from_stack+1} to stack{to_stack+1}"]
                cost = len(new_path) + heuristic(new_stacks, goal_state)
                heapq.heappush(open_set, (cost, new_stacks, new_path))

    return None

initial_state = [['D', 'E', 'H', 'G'], ['A', 'F', 'I'], ['B', 'C']]
goal_state = [['G', 'I'], ['A', 'C', 'D', 'F', 'H'], ['B', 'E']]

solution = a_star_search(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")