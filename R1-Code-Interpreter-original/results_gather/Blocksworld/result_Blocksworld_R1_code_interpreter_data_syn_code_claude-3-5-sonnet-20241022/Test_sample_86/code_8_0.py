def find_solution(initial_stacks, goal_stacks):
    from collections import deque
    
    def get_state_key(stacks):
        return tuple(tuple(stack) for stack in stacks)
    
    def make_move(stacks, from_stack, to_stack):
        if not stacks[from_stack]:
            return None, None
        new_stacks = [list(stack) for stack in stacks]
        block = new_stacks[from_stack].pop()
        new_stacks[to_stack].append(block)
        return new_stacks, block
    
    def is_better_move(current_stacks, from_stack, to_stack, goal_stacks):
        # Check if the move gets a block to its final position
        if not current_stacks[from_stack]:
            return False
        block = current_stacks[from_stack][-1]
        
        # Find where this block should be in goal state
        for stack_idx, stack in enumerate(goal_stacks):
            if stack and block in stack:
                # If we're moving directly to the correct stack
                if to_stack == stack_idx:
                    # Check if it will be in the correct position
                    if len(current_stacks[to_stack]) == stack.index(block):
                        return True
        return False
    
    queue = deque([(initial_stacks, [])])
    visited = {get_state_key(initial_stacks)}
    
    while queue:
        current_stacks, moves = queue.popleft()
        
        # Check if we reached the goal
        if current_stacks == goal_stacks:
            return moves
        
        # First try moves that get blocks to their final positions
        good_moves = []
        other_moves = []
        
        for from_stack in range(len(current_stacks)):
            if not current_stacks[from_stack]:
                continue
            
            for to_stack in range(len(current_stacks)):
                if from_stack != to_stack:
                    new_stacks, block = make_move(current_stacks, from_stack, to_stack)
                    if new_stacks:
                        state_key = get_state_key(new_stacks)
                        if state_key not in visited:
                            move = (new_stacks, moves + [f"Move {block} from {from_stack + 1} to {to_stack + 1}"])
                            if is_better_move(current_stacks, from_stack, to_stack, goal_stacks):
                                good_moves.append(move)
                            else:
                                other_moves.append(move)
        
        # Add good moves first, then other moves
        for new_stacks, new_moves in good_moves + other_moves:
            state_key = get_state_key(new_stacks)
            if state_key not in visited:
                visited.add(state_key)
                queue.append((new_stacks, new_moves))
                
                # Early exit if we found the goal
                if new_stacks == goal_stacks:
                    return new_moves
    
    return None

# Initial and goal states
initial_stacks = [['B', 'C'], ['A'], ['F', 'D', 'E']]
goal_stacks = [['A', 'C'], ['F'], ['B', 'D', 'E']]

# Find and print solution
solution = find_solution(initial_stacks, goal_stacks)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")