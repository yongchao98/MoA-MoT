class State:
    def __init__(self, stacks):
        self.stacks = [list(stack) for stack in stacks]
    
    def copy(self):
        return State([stack.copy() for stack in self.stacks])
    
    def is_valid_move(self, from_stack, to_stack):
        return bool(self.stacks[from_stack])
    
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
    def heuristic(state, goal):
        score = 0
        for i in range(len(state.stacks)):
            for j in range(len(state.stacks[i])):
                if j < len(goal.stacks[i]) and state.stacks[i][j] == goal.stacks[i][j]:
                    score += 1
        return score
    
    def get_possible_moves(state):
        moves = []
        for i in range(len(state.stacks)):
            if state.stacks[i]:
                for j in range(len(state.stacks)):
                    if i != j:
                        new_state = state.copy()
                        if new_state.make_move(i, j):
                            h_score = heuristic(new_state, goal_state)
                            moves.append((h_score, i, j))
        return sorted(moves, reverse=True)
    
    def bfs_with_heuristic(start_state, goal_state):
        queue = [(start_state, [])]
        visited = set()
        
        while queue:
            current_state, path = queue.pop(0)
            state_str = str(current_state)
            
            if current_state == goal_state:
                return path
            
            if state_str in visited:
                continue
                
            visited.add(state_str)
            
            for _, from_stack, to_stack in get_possible_moves(current_state):
                new_state = current_state.copy()
                block = new_state.stacks[from_stack][-1]
                if new_state.make_move(from_stack, to_stack):
                    new_path = path + [(block, from_stack, to_stack)]
                    queue.append((new_state, new_path))
        
        return None

    solution = bfs_with_heuristic(initial_state, goal_state)
    
    if solution:
        moves = []
        for block, from_stack, to_stack in solution:
            moves.append(f"Move {block} from {from_stack+1} to {to_stack+1}")
        print("<<<" + "\n".join(moves) + ">>>")
    else:
        print("No solution found")

initial = State([['A', 'C'], ['E', 'B', 'F', 'D'], []])
goal = State([['B'], ['F'], ['A', 'C', 'D', 'E']])

find_solution(initial, goal)