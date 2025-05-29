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
                for j in range(len(self.stacks)):
                    if i != j:  # can't move to same stack
                        moves.append((i, j, source[-1]))
        return moves
        
    def apply_move(self, move):
        source, dest, _ = move
        if self.stacks[source]:  # if source stack is not empty
            block = self.stacks[source].pop()
            self.stacks[dest].append(block)
            
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
            next_state.apply_move(move)
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
    
    # Format solution
    if solution:
        result = ""
        for block, source, dest in solution:
            result += f"Move {block} from {source} to {dest}\n"
        print(result.strip())

# Define initial and goal states
initial = [
    [],             # stack1
    ['B', 'A', 'C'], # stack2
    ['E', 'D']      # stack3
]

goal = [
    ['B', 'E'],     # stack1
    ['C'],          # stack2
    ['A', 'D']      # stack3
]

solve(initial, goal)