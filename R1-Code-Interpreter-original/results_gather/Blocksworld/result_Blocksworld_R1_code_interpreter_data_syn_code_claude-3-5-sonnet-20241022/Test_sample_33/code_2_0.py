class State:
    def __init__(self, stacks):
        self.stacks = [list(stack) for stack in stacks]
    
    def copy(self):
        return State([stack.copy() for stack in self.stacks])
    
    def is_valid_move(self, from_stack, to_stack, goal_stacks):
        if not self.stacks[from_stack]:
            return False
            
        moving_block = self.stacks[from_stack][-1]
        
        # Find where this block should ultimately go in goal state
        target_stack = None
        for i, stack in enumerate(goal_stacks):
            if moving_block in stack:
                target_stack = i
                break
        
        # If moving to its goal stack, ensure the blocks below are correct
        if to_stack == target_stack:
            goal_stack = goal_stacks[target_stack]
            current_stack = self.stacks[to_stack]
            block_index = goal_stack.index(moving_block)
            
            # Check if blocks below in goal state are in place
            if block_index > 0:
                for i in range(block_index):
                    if i >= len(current_stack) or current_stack[i] != goal_stack[i]:
                        return False
        
        return True
    
    def make_move(self, from_stack, to_stack):
        if self.stacks[from_stack]:
            block = self.stacks[from_stack].pop()
            self.stacks[to_stack].append(block)
            return True
        return False
    
    def __str__(self):
        return str(self.stacks)

def manhattan_distance(current_state, goal_state):
    distance = 0
    block_positions = {}
    
    # Record current positions
    for i, stack in enumerate(current_state.stacks):
        for j, block in enumerate(stack):
            block_positions[block] = (i, j)
    
    # Compare with goal positions
    for i, stack in enumerate(goal_state.stacks):
        for j, block in enumerate(stack):
            if block in block_positions:
                curr_i, curr_j = block_positions[block]
                distance += abs(curr_i - i) + abs(curr_j - j)
    
    return distance

def find_solution(initial_state, goal_state):
    from heapq import heappush, heappop
    
    start = State(initial_state)
    goal = State(goal_state)
    
    # Priority queue with (priority, moves_count, state, moves)
    queue = [(manhattan_distance(start, goal), 0, start, [])]
    visited = {str(start.stacks)}
    
    while queue:
        _, moves_count, current_state, moves = heappop(queue)
        
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
                            new_moves = moves + [move]
                            priority = manhattan_distance(new_state, goal) + len(new_moves)
                            heappush(queue, (priority, len(new_moves), new_state, new_moves))
    
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