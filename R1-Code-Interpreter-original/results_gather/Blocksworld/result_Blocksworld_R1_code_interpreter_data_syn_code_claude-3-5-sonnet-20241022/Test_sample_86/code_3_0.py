class State:
    def __init__(self, stacks):
        self.stacks = stacks
    
    def is_valid_move(self, from_stack, to_stack):
        if (from_stack < 0 or from_stack >= len(self.stacks) or 
            to_stack < 0 or to_stack >= len(self.stacks) or 
            not self.stacks[from_stack]):
            return False
        return True
    
    def make_move(self, from_stack, to_stack):
        if not self.is_valid_move(from_stack, to_stack):
            return None
        new_stacks = [list(stack) for stack in self.stacks]
        block = new_stacks[from_stack].pop()
        new_stacks[to_stack].append(block)
        return State(new_stacks)

def find_solution(initial_state, goal_state):
    from collections import deque
    
    queue = deque([(initial_state, [])])
    visited = {str(initial_state.stacks)}
    
    while queue:
        current_state, moves = queue.popleft()
        
        # Check if we reached the goal
        if current_state.stacks == goal_state.stacks:
            return moves
        
        # Try all possible moves
        for from_stack in range(len(current_state.stacks)):
            if not current_state.stacks[from_stack]:
                continue
            
            # Get the block we're trying to move
            block = current_state.stacks[from_stack][-1]
            
            for to_stack in range(len(current_state.stacks)):
                if from_stack != to_stack:
                    next_state = current_state.make_move(from_stack, to_stack)
                    
                    if next_state and str(next_state.stacks) not in visited:
                        # Create the move description
                        move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                        
                        # Check if this move gets us closer to the goal
                        if next_state.stacks == goal_state.stacks:
                            return moves + [move]
                        
                        visited.add(str(next_state.stacks))
                        queue.append((next_state, moves + [move]))
    
    return None

# Define initial and goal states
initial_stacks = [['B', 'C'], ['A'], ['F', 'D', 'E']]
goal_stacks = [['A', 'C'], ['F'], ['B', 'D', 'E']]

# Create state objects
initial_state = State(initial_stacks)
goal_state = State(goal_stacks)

# Find and print solution
solution = find_solution(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")