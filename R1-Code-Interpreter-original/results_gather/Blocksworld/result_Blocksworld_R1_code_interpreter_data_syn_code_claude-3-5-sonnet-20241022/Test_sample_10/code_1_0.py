def print_moves(initial_state, goal_state):
    # Convert states to lists for easier manipulation
    current = [list(stack) for stack in initial_state]
    goal = [list(stack) for stack in goal_state]
    moves = []
    
    # Helper function to find block location
    def find_block(block, stacks):
        for i, stack in enumerate(stacks):
            if stack and stack[-1] == block:
                return i
        return None

    # Helper function to check if block is accessible (on top)
    def is_accessible(block, stacks):
        for stack in stacks:
            if stack and stack[-1] == block:
                return True
        return False

    # First, build stack1 (G, I)
    # Move G to stack1
    while not (len(current[0]) >= 1 and current[0][-1] == 'G'):
        g_stack = find_block('G', current)
        if is_accessible('G', current):
            moves.append(f"Move G from {g_stack + 1} to 1")
            current[0].append(current[g_stack].pop())

    # Move I to stack1
    while not (len(current[0]) >= 2 and current[0][-1] == 'I'):
        i_stack = find_block('I', current)
        if is_accessible('I', current):
            moves.append(f"Move I from {i_stack + 1} to 1")
            current[0].append(current[i_stack].pop())

    # Build stack2 (A B C D E F H)
    target_stack2 = ['A', 'B', 'C', 'D', 'E', 'F', 'H']
    for block in target_stack2:
        while not (len(current[1]) > 0 and current[1][-1] == block):
            block_stack = find_block(block, current)
            if is_accessible(block, current):
                moves.append(f"Move {block} from {block_stack + 1} to 2")
                current[1].append(current[block_stack].pop())

    # Finally, move J to stack3
    while not (len(current[2]) > 0 and current[2][-1] == 'J'):
        j_stack = find_block('J', current)
        if is_accessible('J', current):
            moves.append(f"Move J from {j_stack + 1} to 3")
            current[2].append(current[j_stack].pop())

    print("<<<" + "\n".join(moves) + ">>>")

# Initial and goal states
initial = [['J'], ['D', 'H', 'F', 'I', 'B', 'A', 'G'], ['E', 'C']]
goal = [['G', 'I'], ['A', 'B', 'C', 'D', 'E', 'F', 'H'], ['J']]

print_moves(initial, goal)