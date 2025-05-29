class State:
    def __init__(self, stacks):
        self.stacks = stacks
        
    def get_top_blocks(self):
        return [stack[-1] if stack else None for stack in self.stacks]
    
    def is_valid_move(self, from_stack, to_stack):
        # Check if source stack has blocks
        if not self.stacks[from_stack]:
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
        
        # Try all possible moves
        for from_stack in range(3):
            for to_stack in range(3):
                if from_stack != to_stack and current_state.is_valid_move(from_stack, to_stack):
                    new_state = current_state.make_move(from_stack, to_stack)
                    if new_state and new_state not in visited:
                        block = current_state.stacks[from_stack][-1]
                        move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                        queue.append((new_state, moves + [move]))
                        visited.add(new_state)
    
    return None

# Initial and goal states
initial = State([['A'], ['E', 'C'], ['D', 'B']])
goal = State([['B', 'D'], ['E'], ['A', 'C']])

# Find solution
solution = find_solution(initial, goal)
if solution:
    print("<<<")
    print("\n".join(solution))
    print(">>>")
else:
    print("No solution found")