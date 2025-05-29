class State:
    def __init__(self, stacks):
        self.stacks = [list(stack) for stack in stacks]
    
    def copy(self):
        return State([stack.copy() for stack in self.stacks])
    
    def is_valid_move(self, from_stack, to_stack, goal_stacks):
        # Check if source stack has blocks
        if not self.stacks[from_stack]:
            return False
        
        # Get the block to be moved
        block = self.stacks[from_stack][-1]
        
        # Check if this move would place the block in a position that makes the goal impossible
        # Find where this block should be in the goal state
        goal_stack_idx = None
        for i, stack in enumerate(goal_stacks):
            if block in stack:
                goal_stack_idx = i
                break
        
        # If we're moving to a stack other than the goal stack, make sure it's a temporary move
        if to_stack != goal_stack_idx:
            # Allow temporary moves
            return True
        
        # If we're moving to the goal stack, check if it's the right time
        if to_stack == goal_stack_idx:
            goal_stack = goal_stacks[to_stack]
            block_idx = goal_stack.index(block)
            # Check if all blocks below this one in goal state are already in place
            current_stack = self.stacks[to_stack]
            if block_idx > 0:
                for i in range(block_idx):
                    if i >= len(current_stack) or current_stack[i] != goal_stack[i]:
                        return False
        
        return True
    
    def make_move(self, from_stack, to_stack):
        if self.stacks[from_stack]:
            block = self.stacks[from_stack].pop()
            self.stacks[to_stack].append(block)
            return True
        return False

def find_solution(initial_state, goal_state):
    from collections import deque
    
    start = State(initial_state)
    goal = State(goal_state)
    queue = deque([(start, [])])
    visited = {str(start.stacks)}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if current_state.stacks == goal.stacks:
            return moves
        
        for from_stack in range(len(current_state.stacks)):
            for to_stack in range(len(current_state.stacks)):
                if from_stack != to_stack and current_state.stacks[from_stack]:
                    new_state = current_state.copy()
                    if new_state.is_valid_move(from_stack, to_stack, goal.stacks) and new_state.make_move(from_stack, to_stack):
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