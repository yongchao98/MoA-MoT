class State:
    def __init__(self, stacks):
        self.stacks = [list(stack) for stack in stacks]
    
    def get_top_block(self, stack_idx):
        if 0 <= stack_idx < len(self.stacks) and self.stacks[stack_idx]:
            return self.stacks[stack_idx][-1]
        return None

def find_solution(initial_stacks, goal_stacks):
    def get_state_key(stacks):
        return tuple(tuple(stack) for stack in stacks)
    
    def make_move(stacks, from_stack, to_stack):
        new_stacks = [list(stack) for stack in stacks]
        block = new_stacks[from_stack].pop()
        new_stacks[to_stack].append(block)
        return new_stacks, block
    
    from collections import deque
    
    # Initialize queue with (stacks, moves)
    queue = deque([(initial_stacks, [])])
    visited = {get_state_key(initial_stacks)}
    
    while queue:
        current_stacks, moves = queue.popleft()
        
        # Check if we reached the goal
        if current_stacks == goal_stacks:
            return moves
        
        # Try all possible moves
        for from_stack in range(len(current_stacks)):
            # Skip if source stack is empty
            if not current_stacks[from_stack]:
                continue
            
            # Try moving to each other stack
            for to_stack in range(len(current_stacks)):
                if from_stack != to_stack:
                    # Make the move
                    try:
                        new_stacks, block = make_move(current_stacks, from_stack, to_stack)
                        state_key = get_state_key(new_stacks)
                        
                        if state_key not in visited:
                            visited.add(state_key)
                            new_moves = moves + [f"Move {block} from {from_stack + 1} to {to_stack + 1}"]
                            queue.append((new_stacks, new_moves))
                            
                            # Early exit if we found the goal
                            if new_stacks == goal_stacks:
                                return new_moves
                    except:
                        continue
    
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