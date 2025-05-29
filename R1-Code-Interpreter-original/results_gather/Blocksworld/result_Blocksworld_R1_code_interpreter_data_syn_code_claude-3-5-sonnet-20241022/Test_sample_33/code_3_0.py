class State:
    def __init__(self, stacks):
        self.stacks = [list(stack) for stack in stacks]
    
    def copy(self):
        return State([stack.copy() for stack in self.stacks])
    
    def is_valid_move(self, from_stack, to_stack):
        return len(self.stacks[from_stack]) > 0
    
    def make_move(self, from_stack, to_stack):
        if self.is_valid_move(from_stack, to_stack):
            block = self.stacks[from_stack].pop()
            self.stacks[to_stack].append(block)
            return True
        return False
    
    def __str__(self):
        return str(self.stacks)

def find_solution(initial_state, goal_state):
    from collections import deque
    
    def is_goal_reached(current, goal):
        return all(current.stacks[i] == goal.stacks[i] for i in range(len(current.stacks)))
    
    start = State(initial_state)
    goal = State(goal_state)
    queue = deque([(start, [])])
    visited = {str(start.stacks)}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if is_goal_reached(current_state, goal):
            return moves
        
        for from_stack in range(len(current_state.stacks)):
            if not current_state.stacks[from_stack]:  # Skip empty stacks
                continue
                
            for to_stack in range(len(current_state.stacks)):
                if from_stack != to_stack:
                    new_state = current_state.copy()
                    if new_state.make_move(from_stack, to_stack):
                        move = f"Move {current_state.stacks[from_stack][-1]} from {from_stack + 1} to {to_stack + 1}"
                        state_str = str(new_state.stacks)
                        
                        if state_str not in visited:
                            visited.add(state_str)
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