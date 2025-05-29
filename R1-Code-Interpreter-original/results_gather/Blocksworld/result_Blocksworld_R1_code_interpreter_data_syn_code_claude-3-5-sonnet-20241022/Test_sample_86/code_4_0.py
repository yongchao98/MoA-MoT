class State:
    def __init__(self, stacks):
        self.stacks = [list(stack) for stack in stacks]
    
    def is_valid_move(self, from_stack, to_stack):
        return (0 <= from_stack < len(self.stacks) and 
                0 <= to_stack < len(self.stacks) and 
                len(self.stacks[from_stack]) > 0)
    
    def make_move(self, from_stack, to_stack):
        if not self.is_valid_move(from_stack, to_stack):
            return None
        new_state = State(self.stacks)
        block = new_state.stacks[from_stack].pop()
        new_state.stacks[to_stack].append(block)
        return new_state
    
    def evaluate(self, goal_state):
        # Count how many blocks are in their correct final position
        score = 0
        for i in range(len(self.stacks)):
            current = self.stacks[i]
            goal = goal_state.stacks[i]
            for j in range(min(len(current), len(goal))):
                if current[j] == goal[j]:
                    score += 1
        return score

def find_solution(initial_state, goal_state):
    from collections import deque
    
    best_solution = None
    best_length = float('inf')
    
    def dfs(current_state, moves, visited):
        nonlocal best_solution, best_length
        
        # If we've found a solution
        if current_state.stacks == goal_state.stacks:
            if len(moves) < best_length:
                best_solution = moves.copy()
                best_length = len(moves)
            return
        
        # If we've made too many moves, backtrack
        if len(moves) >= best_length:
            return
        
        # Get current evaluation score
        current_score = current_state.evaluate(goal_state)
        
        # Try all possible moves
        possible_moves = []
        for from_stack in range(len(current_state.stacks)):
            if not current_state.stacks[from_stack]:
                continue
            for to_stack in range(len(current_state.stacks)):
                if from_stack != to_stack:
                    next_state = current_state.make_move(from_stack, to_stack)
                    if next_state:
                        next_score = next_state.evaluate(goal_state)
                        state_str = str(next_state.stacks)
                        if state_str not in visited:
                            possible_moves.append((next_score, from_stack, to_stack, next_state))
        
        # Sort moves by evaluation score (higher is better)
        possible_moves.sort(reverse=True)
        
        # Try each move
        for _, from_stack, to_stack, next_state in possible_moves:
            block = current_state.stacks[from_stack][-1]
            move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
            
            state_str = str(next_state.stacks)
            if state_str not in visited:
                visited.add(state_str)
                moves.append(move)
                dfs(next_state, moves, visited)
                moves.pop()
                visited.remove(state_str)
    
    # Start the search
    visited = {str(initial_state.stacks)}
    dfs(initial_state, [], visited)
    
    return best_solution

# Initial and goal states
initial_stacks = [['B', 'C'], ['A'], ['F', 'D', 'E']]
goal_stacks = [['A', 'C'], ['F'], ['B', 'D', 'E']]

# Create state objects
initial_state = State(initial_stacks)
goal_state = State(goal_stacks)

# Find and print solution
solution = find_solution(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")