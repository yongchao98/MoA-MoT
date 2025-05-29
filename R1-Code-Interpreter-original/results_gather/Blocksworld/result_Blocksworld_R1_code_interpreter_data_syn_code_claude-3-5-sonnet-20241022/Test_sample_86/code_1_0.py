def find_direct_solution(initial_stacks, goal_stacks):
    # We know the exact sequence needed:
    # 1. Move C from stack1 to stack2
    # 2. Move B from stack1 to stack3
    # 3. Move A from stack2 to stack1
    # 4. Move C from stack2 to stack1
    # 5. Move F from stack3 to stack2
    
    solution = []
    current_stacks = [list(stack) for stack in initial_stacks]
    
    # Helper function to find the top block of a stack
    def get_top_block(stack_idx):
        if current_stacks[stack_idx]:
            return current_stacks[stack_idx][-1]
        return None
    
    # Helper function to make a move
    def make_move(from_stack, to_stack):
        if current_stacks[from_stack]:
            block = current_stacks[from_stack].pop()
            current_stacks[to_stack].append(block)
            return f"Move {block} from {from_stack + 1} to {to_stack + 1}"
        return None
    
    # Step 1: Move C from stack1 to stack2
    if get_top_block(0) == 'C':
        move = make_move(0, 1)
        if move:
            solution.append(move)
    
    # Step 2: Move B from stack1 to stack3
    if get_top_block(0) == 'B':
        move = make_move(0, 2)
        if move:
            solution.append(move)
    
    # Step 3: Move A from stack2 to stack1
    if get_top_block(1) == 'A':
        move = make_move(1, 0)
        if move:
            solution.append(move)
    
    # Step 4: Move C from stack2 to stack1
    if get_top_block(1) == 'C':
        move = make_move(1, 0)
        if move:
            solution.append(move)
    
    # Step 5: Move F from stack3 to stack2
    if get_top_block(2) == 'F':
        move = make_move(2, 1)
        if move:
            solution.append(move)
    
    # Verify that we reached the goal state
    if current_stacks == goal_stacks:
        return solution
    return None

# Initial and goal states
initial_stacks = [['B', 'C'], ['A'], ['F', 'D', 'E']]
goal_stacks = [['A', 'C'], ['F'], ['B', 'D', 'E']]

# Find and print solution
solution = find_direct_solution(initial_stacks, goal_stacks)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")