class State:
    def __init__(self, stacks):
        self.stacks = [list(stack) for stack in stacks]
    
    def copy(self):
        return State([stack.copy() for stack in self.stacks])
    
    def is_valid_move(self, from_stack, to_stack):
        # Check if source stack has blocks
        if not self.stacks[from_stack]:
            return False
        return True
    
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
    
    start = State(initial_state)
    goal = State(goal_state)
    queue = deque([(start, [])])
    visited = {start}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if current_state == goal:
            return moves
        
        # Try all possible moves
        for from_stack in range(len(current_state.stacks)):
            for to_stack in range(len(current_state.stacks)):
                if from_stack != to_stack and current_state.stacks[from_stack]:
                    new_state = current_state.copy()
                    if new_state.make_move(from_stack, to_stack):
                        move = f"Move {current_state.stacks[from_stack][-1]} from {from_stack + 1} to {to_stack + 1}"
                        if new_state not in visited:
                            visited.add(new_state)
                            queue.append((new_state, moves + [move]))
    
    return None

# Initial and goal states
initial = [['A'], ['E', 'B', 'D', 'F'], ['G', 'C'], []]
goal = [['A', 'D', 'E', 'F'], ['G'], ['B', 'C'], []]

# Find solution
solution = find_solution(initial, goal)
if solution:
    print('\n'.join(solution))
else:
    print("No solution found")