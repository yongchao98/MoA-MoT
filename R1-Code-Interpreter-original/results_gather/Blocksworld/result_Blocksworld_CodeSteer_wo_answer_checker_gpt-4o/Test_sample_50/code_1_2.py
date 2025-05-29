def is_goal_state(stacks, goal_state):
    return stacks == goal_state

def move_block(stacks, from_stack, to_stack):
    new_stacks = [list(stack) for stack in stacks]
    block = new_stacks[from_stack].pop()
    new_stacks[to_stack].append(block)
    return new_stacks

def solve_blocks(stacks, goal_state, path):
    if is_goal_state(stacks, goal_state):
        return path

    for from_stack in range(len(stacks)):
        if not stacks[from_stack]:
            continue
        for to_stack in range(len(stacks)):
            if from_stack == to_stack:
                continue
            new_stacks = move_block(stacks, from_stack, to_stack)
            new_path = path + [f"Move {stacks[from_stack][-1]} from stack{from_stack+1} to stack{to_stack+1}"]
            result = solve_blocks(new_stacks, goal_state, new_path)
            if result:
                return result
    return None

initial_state = [['D', 'E', 'H', 'G'], ['A', 'F', 'I'], ['B', 'C']]
goal_state = [['G', 'I'], ['A', 'C', 'D', 'F', 'H'], ['B', 'E']]

solution = solve_blocks(initial_state, goal_state, [])
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")