class State:
    def __init__(self, stacks):
        self.stacks = stacks
        
    def get_top_blocks(self):
        return [stack[-1] if stack else None for stack in self.stacks]
    
    def is_valid_move(self, from_stack, to_stack):
        if from_stack < 0 or from_stack >= len(self.stacks) or not self.stacks[from_stack]:
            return False
        return True
    
    def make_move(self, from_stack, to_stack):
        if not self.is_valid_move(from_stack, to_stack):
            return None
        new_stacks = [list(stack) for stack in self.stacks]
        block = new_stacks[from_stack].pop()
        new_stacks[to_stack].append(block)
        return State(new_stacks)
    
    def __eq__(self, other):
        return self.stacks == other.stacks
    
    def __hash__(self):
        return hash(tuple(tuple(stack) for stack in self.stacks))

def find_solution(initial_state, goal_state):
    from collections import deque, defaultdict
    
    queue = deque([(initial_state, [])])
    visited = {initial_state}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if current_state == goal_state:
            return moves
        
        top_blocks = current_state.get_top_blocks()
        
        for from_stack in range(len(current_state.stacks)):
            if not current_state.stacks[from_stack]:
                continue
                
            for to_stack in range(len(current_state.stacks)):
                if from_stack != to_stack:
                    next_state = current_state.make_move(from_stack, to_stack)
                    if next_state and next_state not in visited:
                        move = f"Move {top_blocks[from_stack]} from {from_stack + 1} to {to_stack + 1}"
                        visited.add(next_state)
                        queue.append((next_state, moves + [move]))
    
    return None

# Initial and goal states
initial_stacks = [['B', 'C'], ['A'], ['F', 'D', 'E']]
goal_stacks = [['A', 'C'], ['F'], ['B', 'D', 'E']]

initial_state = State(initial_stacks)
goal_state = State(goal_stacks)

# Find solution
solution = find_solution(initial_state, goal_state)

if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")