class State:
    def __init__(self, stacks):
        self.stacks = stacks
        
    def is_valid_move(self, from_stack, to_stack):
        # Check if source stack is not empty
        if not self.stacks[from_stack]:
            return False
        return True
    
    def make_move(self, from_stack, to_stack):
        if self.is_valid_move(from_stack, to_stack):
            block = self.stacks[from_stack][-1]  # Get top block
            new_stacks = [stack[:] for stack in self.stacks]  # Deep copy
            new_stacks[from_stack].pop()  # Remove from source
            new_stacks[to_stack].append(block)  # Add to destination
            return State(new_stacks), f"Move {block} from {from_stack+1} to {to_stack+1}"
        return None, None

def is_goal_state(state):
    goal = [['A', 'E'], ['B'], ['C', 'D'], []]
    return state.stacks == goal

def solve_blocksworld():
    initial_stacks = [['E', 'D'], ['A', 'B'], ['C'], []]
    initial_state = State(initial_stacks)
    visited = set()
    queue = [(initial_state, [])]
    
    while queue:
        current_state, moves = queue.pop(0)
        state_tuple = tuple(tuple(stack) for stack in current_state.stacks)
        
        if state_tuple in visited:
            continue
            
        visited.add(state_tuple)
        
        if is_goal_state(current_state):
            return moves
            
        for from_stack in range(4):
            for to_stack in range(4):
                if from_stack != to_stack:
                    new_state, move = current_state.make_move(from_stack, to_stack)
                    if new_state:
                        new_moves = moves + [move]
                        queue.append((new_state, new_moves))
    
    return None

# Solve and print the solution
solution = solve_blocksworld()
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")