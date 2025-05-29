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
    
    def __str__(self):
        return str(self.stacks)

def find_solution(initial_state, goal_state):
    def get_possible_moves(state):
        moves = []
        for i in range(len(state.stacks)):
            if state.stacks[i]:  # if stack has blocks
                for j in range(len(state.stacks)):
                    if i != j:  # don't move to same stack
                        moves.append((i, j))
        return moves
    
    def dfs(current_state, goal_state, path, visited, max_depth=20):
        if len(path) > max_depth:
            return None
        if current_state == goal_state:
            return path
        
        state_str = str(current_state)
        if state_str in visited and len(path) >= visited[state_str]:
            return None
        visited[state_str] = len(path)
        
        for from_stack, to_stack in get_possible_moves(current_state):
            if current_state.stacks[from_stack]:  # Check if source stack has blocks
                new_state = current_state.copy()
                block = new_state.stacks[from_stack][-1]  # Get top block
                if new_state.make_move(from_stack, to_stack):
                    new_path = path + [(block, from_stack, to_stack)]
                    result = dfs(new_state, goal_state, new_path, visited)
                    if result is not None:
                        return result
        return None

    # Initialize states
    visited = {}
    path = []
    solution = dfs(initial_state, goal_state, path, visited)
    
    if solution:
        moves = []
        for block, from_stack, to_stack in solution:
            moves.append(f"Move {block} from {from_stack+1} to {to_stack+1}")
        print("<<<" + "\n".join(moves) + ">>>")
    else:
        print("No solution found")

# Initial and goal states
initial = State([['A', 'C'], ['E', 'B', 'F', 'D'], []])
goal = State([['B'], ['F'], ['A', 'C', 'D', 'E']])

find_solution(initial, goal)