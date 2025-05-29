class State:
    def __init__(self, stacks):
        self.stacks = [list(stack) for stack in stacks]
    
    def copy(self):
        return State([stack.copy() for stack in self.stacks])
    
    def is_valid_move(self, from_stack, to_stack):
        return (from_stack >= 0 and from_stack < len(self.stacks) and 
                to_stack >= 0 and to_stack < len(self.stacks) and 
                len(self.stacks[from_stack]) > 0)
    
    def make_move(self, from_stack, to_stack):
        if not self.is_valid_move(from_stack, to_stack):
            return None
        new_state = self.copy()
        block = new_state.stacks[from_stack].pop()
        new_state.stacks[to_stack].append(block)
        return new_state
    
    def __eq__(self, other):
        return self.stacks == other.stacks
    
    def __hash__(self):
        return hash(tuple(tuple(stack) for stack in self.stacks))

def find_solution(initial_state, goal_state):
    from collections import deque
    
    queue = deque([(initial_state, [])])
    visited = {initial_state}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if current_state == goal_state:
            return moves
        
        for from_stack in range(len(current_state.stacks)):
            for to_stack in range(len(current_state.stacks)):
                if from_stack != to_stack and current_state.is_valid_move(from_stack, to_stack):
                    next_state = current_state.make_move(from_stack, to_stack)
                    if next_state not in visited:
                        block = current_state.stacks[from_stack][-1]
                        move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                        queue.append((next_state, moves + [move]))
                        visited.add(next_state)
    
    return None

# Initial and goal states
initial = State([['C', 'A', 'F'], ['B'], ['D', 'E']])
goal = State([['B'], ['A', 'C', 'D'], ['E', 'F']])

# Find solution
solution = find_solution(initial, goal)
if solution:
    print('\n'.join(solution))