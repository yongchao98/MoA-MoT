class State:
    def __init__(self, stacks):
        self.stacks = stacks
    
    def is_valid_move(self, from_stack, to_stack):
        # Check bounds and if source stack has blocks
        if (from_stack < 0 or from_stack >= len(self.stacks) or 
            to_stack < 0 or to_stack >= len(self.stacks) or 
            not self.stacks[from_stack]):
            return False
        # Can only move top block
        return True
    
    def make_move(self, from_stack, to_stack):
        if not self.is_valid_move(from_stack, to_stack):
            return None
        new_stacks = [list(stack) for stack in self.stacks]
        block = new_stacks[from_stack].pop()
        new_stacks[to_stack].append(block)
        return State(new_stacks)
    
    def __eq__(self, other):
        return str(self.stacks) == str(other.stacks)
    
    def __hash__(self):
        return hash(str(self.stacks))
    
    def manhattan_distance(self, goal_state):
        # Heuristic function to guide the search
        distance = 0
        for i, stack in enumerate(self.stacks):
            for j, block in enumerate(stack):
                # Find where this block should be in goal state
                for gi, gstack in enumerate(goal_state.stacks):
                    if block in gstack:
                        distance += abs(i - gi)
        return distance

def find_solution(initial_state, goal_state):
    from heapq import heappush, heappop
    
    # Priority queue with (priority, moves_count, state, moves)
    queue = [(0, 0, initial_state, [])]
    visited = {str(initial_state.stacks)}
    moves_count = 0
    
    while queue:
        _, _, current_state, moves = heappop(queue)
        
        if current_state.stacks == goal_state.stacks:
            return moves
        
        moves_count += 1
        for from_stack in range(len(current_state.stacks)):
            if not current_state.stacks[from_stack]:
                continue
            
            for to_stack in range(len(current_state.stacks)):
                if from_stack != to_stack:
                    next_state = current_state.make_move(from_stack, to_stack)
                    if next_state and str(next_state.stacks) not in visited:
                        move = f"Move {current_state.stacks[from_stack][-1]} from {from_stack + 1} to {to_stack + 1}"
                        new_moves = moves + [move]
                        visited.add(str(next_state.stacks))
                        # Use manhattan distance as heuristic
                        priority = len(new_moves) + next_state.manhattan_distance(goal_state)
                        heappush(queue, (priority, moves_count, next_state, new_moves))
    
    return None

# Initial and goal states
initial_stacks = [['B', 'C'], ['A'], ['F', 'D', 'E']]
goal_stacks = [['A', 'C'], ['F'], ['B', 'D', 'E']]

initial_state = State(initial_stacks)
goal_state = State(goal_stacks)

# Find and print solution
solution = find_solution(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")