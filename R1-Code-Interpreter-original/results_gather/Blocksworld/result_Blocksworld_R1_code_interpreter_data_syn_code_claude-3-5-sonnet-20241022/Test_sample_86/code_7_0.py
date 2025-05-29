def get_block_position(stacks, block):
    for i, stack in enumerate(stacks):
        if block in stack:
            return i
    return -1

def is_move_towards_goal(current_stacks, from_stack, to_stack, block, goal_stacks):
    current_pos = from_stack
    goal_pos = get_block_position(goal_stacks, block)
    
    # If moving to the goal position directly
    if to_stack == goal_pos:
        return True
    
    # If current position is wrong and we're moving away from it
    if current_pos != goal_pos:
        return True
        
    return False

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
    
    queue = deque([(initial_stacks, [])])
    visited = {get_state_key(initial_stacks)}
    
    # Keep track of minimum solution length found
    min_solution = None
    min_length = float('inf')
    
    while queue:
        current_stacks, moves = queue.popleft()
        
        # Skip if we already found a shorter solution
        if len(moves) >= min_length:
            continue
        
        # Check if we reached the goal
        if current_stacks == goal_stacks:
            if len(moves) < min_length:
                min_solution = moves
                min_length = len(moves)
            continue
        
        # Try all possible moves
        possible_moves = []
        for from_stack in range(len(current_stacks)):
            if not current_stacks[from_stack]:
                continue
            
            block = current_stacks[from_stack][-1]
            for to_stack in range(len(current_stacks)):
                if from_stack != to_stack:
                    new_stacks, moved_block = make_move(current_stacks, from_stack, to_stack)
                    if new_stacks:
                        state_key = get_state_key(new_stacks)
                        if state_key not in visited:
                            # Prioritize moves that get blocks to their goal positions
                            priority = is_move_towards_goal(current_stacks, from_stack, to_stack, block, goal_stacks)
                            possible_moves.append((priority, from_stack, to_stack, new_stacks, moved_block))
        
        # Sort moves by priority (True comes before False)
        possible_moves.sort(reverse=True)
        
        # Try each move
        for _, from_stack, to_stack, new_stacks, block in possible_moves:
            state_key = get_state_key(new_stacks)
            if state_key not in visited:
                visited.add(state_key)
                new_moves = moves + [f"Move {block} from {from_stack + 1} to {to_stack + 1}"]
                queue.append((new_stacks, new_moves))
    
    return min_solution

# Initial and goal states
initial_stacks = [['B', 'C'], ['A'], ['F', 'D', 'E']]
goal_stacks = [['A', 'C'], ['F'], ['B', 'D', 'E']]

# Find and print solution
solution = find_solution(initial_stacks, goal_stacks)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")