def is_goal_state(stacks, goal_state):
    return stacks == goal_state

def move_block(stacks, from_stack, to_stack):
    new_stacks = [list(stack) for stack in stacks]
    block = new_stacks[from_stack].pop()
    new_stacks[to_stack].append(block)
    return new_stacks

def stacks_to_tuple(stacks):
    return tuple(tuple(stack) for stack in stacks)

def iterative_deepening_search(initial_state, goal_state):
    def dls(stacks, path, depth):
        if is_goal_state(stacks, goal_state):
            return path
        if depth == 0:
            return None
        visited.add(stacks_to_tuple(stacks))
        for from_stack in range(len(stacks)):
            if not stacks[from_stack]:
                continue
            for to_stack in range(len(stacks)):
                if from_stack == to_stack:
                    continue
                new_stacks = move_block(stacks, from_stack, to_stack)
                new_stacks_tuple = stacks_to_tuple(new_stacks)
                if new_stacks_tuple not in visited:
                    new_path = path + [f"Move {stacks[from_stack][-1]} from stack{from_stack+1} to stack{to_stack+1}"]
                    result = dls(new_stacks, new_path, depth - 1)
                    if result:
                        return result
        visited.remove(stacks_to_tuple(stacks))
        return None

    depth = 0
    while True:
        visited = set()
        result = dls(initial_state, [], depth)
        if result:
            return result
        depth += 1

initial_state = [['D', 'E', 'H', 'G'], ['A', 'F', 'I'], ['B', 'C']]
goal_state = [['G', 'I'], ['A', 'C', 'D', 'F', 'H'], ['B', 'E']]

solution = iterative_deepening_search(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")