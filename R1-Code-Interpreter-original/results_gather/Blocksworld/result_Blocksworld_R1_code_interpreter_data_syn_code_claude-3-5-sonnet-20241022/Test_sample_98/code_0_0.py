class State:
    def __init__(self, stacks):
        self.stacks = stacks
        
    def copy(self):
        return State([stack[:] for stack in self.stacks])
        
    def is_goal(self, goal):
        return all(self.stacks[i] == goal.stacks[i] for i in range(len(self.stacks)))
        
    def get_valid_moves(self):
        moves = []
        for i, source in enumerate(self.stacks):
            if source:  # if stack is not empty
                for j, dest in enumerate(self.stacks):
                    if i != j:  # can't move to same stack
                        moves.append((i, j, source[-1]))
        return moves
        
    def apply_move(self, move):
        source, dest, block = move
        if self.stacks[source] and self.stacks[source][-1] == block:
            self.stacks[dest].append(self.stacks[source].pop())
            return True
        return False

def solve(initial, goal):
    def dfs(current, path, visited):
        if current.is_goal(goal):
            return path
            
        state_str = str(current.stacks)
        if state_str in visited:
            return None
        visited.add(state_str)
        
        for move in current.get_valid_moves():
            next_state = current.copy()
            if next_state.apply_move(move):
                new_path = path + [(move[2], move[0]+1, move[1]+1)]
                result = dfs(next_state, new_path, visited)
                if result is not None:
                    return result
        return None

    # Initial and goal states
    initial_state = State(initial)
    goal_state = State(goal)
    
    # Find solution
    solution = dfs(initial_state, [], set())
    
    if solution:
        # Format solution
        formatted = []
        for block, source, dest in solution:
            formatted.append(f"Move {block} from {source} to {dest}")
        print("<<<" + "\n".join(formatted) + ">>>")
    else:
        print("No solution found")

# Define initial and goal states
initial = [['A', 'E'], ['D'], ['C', 'B']]
goal = [['A', 'B', 'E'], ['C'], ['D']]

solve(initial, goal)