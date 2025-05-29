class State:
    def __init__(self, stacks):
        self.stacks = [list(stack) for stack in stacks]
    
    def copy(self):
        return State([stack.copy() for stack in self.stacks])
    
    def is_valid_move(self, from_stack, to_stack):
        return (0 <= from_stack < len(self.stacks) and 
                0 <= to_stack < len(self.stacks) and 
                len(self.stacks[from_stack]) > 0)
    
    def make_move(self, from_stack, to_stack):
        if self.is_valid_move(from_stack, to_stack):
            block = self.stacks[from_stack].pop()
            self.stacks[to_stack].append(block)
            return True
        return False
    
    def __eq__(self, other):
        return all(self.stacks[i] == other.stacks[i] for i in range(len(self.stacks)))
    
    def __hash__(self):
        return hash(tuple(tuple(stack) for stack in self.stacks))

def find_solution(initial_state, goal_state):
    from collections import deque
    
    queue = deque([(initial_state, [])])
    visited = set()
    
    while queue:
        current_state, moves = queue.popleft()
        state_hash = hash(current_state)
        
        if current_state == goal_state:
            return moves
        
        if state_hash in visited:
            continue
            
        visited.add(state_hash)
        
        for from_stack in range(len(current_state.stacks)):
            for to_stack in range(len(current_state.stacks)):
                if from_stack != to_stack and current_state.is_valid_move(from_stack, to_stack):
                    new_state = current_state.copy()
                    block = new_state.stacks[from_stack][-1]
                    if new_state.make_move(from_stack, to_stack):
                        new_moves = moves + [f"Move {block} from {from_stack+1} to {to_stack+1}"]
                        queue.append((new_state, new_moves))

# Initial and goal states
initial = State([['G', 'B', 'F', 'E', 'D', 'A'], [], ['C']])
goal = State([['C', 'G'], ['A', 'D', 'E'], ['B', 'F']])

# Find and print solution
solution = find_solution(initial, goal)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")