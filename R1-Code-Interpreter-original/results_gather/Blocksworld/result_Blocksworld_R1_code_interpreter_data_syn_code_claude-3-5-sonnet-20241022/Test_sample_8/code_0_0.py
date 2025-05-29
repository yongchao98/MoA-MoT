from collections import deque
import copy

class State:
    def __init__(self, stacks):
        self.stacks = stacks
        
    def __str__(self):
        return str(self.stacks)
    
    def __eq__(self, other):
        return str(self) == str(other)
    
    def __hash__(self):
        return hash(str(self))

def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(state, from_stack, to_stack):
    if from_stack < 0 or from_stack >= len(state.stacks):
        return False
    if to_stack < 0 or to_stack >= len(state.stacks):
        return False
    if not state.stacks[from_stack]:
        return False
    return True

def make_move(state, from_stack, to_stack):
    if not is_valid_move(state, from_stack, to_stack):
        return None
    new_stacks = copy.deepcopy(state.stacks)
    block = new_stacks[from_stack].pop()
    new_stacks[to_stack].append(block)
    return State(new_stacks)

def find_solution():
    initial_stacks = [['D'], ['G', 'C', 'A', 'F'], ['B', 'I', 'E', 'H']]
    goal_stacks = [['A', 'F', 'H'], ['C', 'G', 'I'], ['B', 'D', 'E']]
    
    initial_state = State(initial_stacks)
    goal_state = State(goal_stacks)
    
    queue = deque([(initial_state, [])])
    visited = {initial_state}
    
    while queue:
        current_state, moves = queue.popleft()
        
        if current_state == goal_state:
            return moves
        
        for from_stack in range(3):
            for to_stack in range(3):
                if from_stack != to_stack:
                    if current_state.stacks[from_stack]:
                        new_state = make_move(current_state, from_stack, to_stack)
                        if new_state and new_state not in visited:
                            block = current_state.stacks[from_stack][-1]
                            new_moves = moves + [f"Move {block} from {from_stack + 1} to {to_stack + 1}"]
                            queue.append((new_state, new_moves))
                            visited.add(new_state)
    return None

solution = find_solution()
if solution:
    print("<<<")
    print("\n".join(solution))
    print(">>>")
else:
    print("No solution found")