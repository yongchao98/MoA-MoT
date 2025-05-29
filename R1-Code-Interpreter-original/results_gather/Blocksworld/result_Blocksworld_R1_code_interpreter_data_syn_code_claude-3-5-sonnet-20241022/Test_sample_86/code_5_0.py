class State:
    def __init__(self, stacks):
        self.stacks = [list(stack) for stack in stacks]
    
    def is_valid_move(self, from_stack, to_stack):
        return (0 <= from_stack < len(self.stacks) and 
                0 <= to_stack < len(self.stacks) and 
                len(self.stacks[from_stack]) > 0 and
                from_stack != to_stack)
    
    def make_move(self, from_stack, to_stack):
        if not self.is_valid_move(from_stack, to_stack):
            return None
        new_state = State(self.stacks)
        block = new_state.stacks[from_stack].pop()
        new_state.stacks[to_stack].append(block)
        return new_state, block

def find_solution(initial_state, goal_state):
    from collections import deque
    
    queue = deque([(initial_state, [])])
    visited = {str(initial_state.stacks)}
    max_moves = 10  # Limit the number of moves to ensure quick solution
    
    while queue:
        current_state, moves = queue.popleft()
        
        if len(moves) >= max_moves:
            continue
            
        if current_state.stacks == goal_state.stacks:
            return moves
            
        for from_stack in range(len(current_state.stacks)):
            if not current_state.stacks[from_stack]:
                continue
                
            for to_stack in range(len(current_state.stacks)):
                if from_stack == to_stack:
                    continue
                    
                result = current_state.make_move(from_stack, to_stack)
                if result:
                    next_state, block = result
                    state_str = str(next_state.stacks)
                    
                    if state_str not in visited:
                        move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                        visited.add(state_str)
                        queue.append((next_state, moves + [move]))
                        
                        # Early exit if we found the goal
                        if next_state.stacks == goal_state.stacks:
                            return moves + [move]
    
    return None

# Initial and goal states
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