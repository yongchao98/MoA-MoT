class State:
    def __init__(self, stacks):
        self.stacks = [list(stack) for stack in stacks]

    def is_valid_move(self, from_stack, to_stack):
        return len(self.stacks[from_stack-1]) > 0

    def move(self, from_stack, to_stack):
        if not self.is_valid_move(from_stack, to_stack):
            return None
        new_stacks = [list(stack) for stack in self.stacks]
        block = new_stacks[from_stack-1].pop()
        new_stacks[to_stack-1].append(block)
        return State(new_stacks)

    def __str__(self):
        return str(self.stacks)

    def __eq__(self, other):
        return str(self) == str(other)

def find_solution(initial_state, goal_state):
    from collections import deque
    
    visited = set()
    queue = deque([(initial_state, [])])
    
    while queue:
        current_state, moves = queue.popleft()
        state_str = str(current_state)
        
        if current_state == goal_state:
            return moves
            
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        for from_stack in range(1, 4):
            if len(current_state.stacks[from_stack-1]) == 0:
                continue
            block = current_state.stacks[from_stack-1][-1]
            for to_stack in range(1, 4):
                if from_stack != to_stack:
                    next_state = current_state.move(from_stack, to_stack)
                    if next_state:
                        move = f"Move {block} from {from_stack} to {to_stack}"
                        queue.append((next_state, moves + [move]))

# Initial and goal states
initial = State([['A', 'C'], ['E', 'D', 'B', 'H', 'F'], ['G']])
goal = State([['A', 'G', 'H'], ['F'], ['B', 'C', 'D', 'E']])

# Find solution
solution = find_solution(initial, goal)
if solution:
    print('\n'.join(solution))