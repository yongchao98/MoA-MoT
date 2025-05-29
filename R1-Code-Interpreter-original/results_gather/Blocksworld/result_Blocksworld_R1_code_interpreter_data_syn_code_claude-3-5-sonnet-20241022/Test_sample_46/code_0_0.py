from collections import deque
import copy

class BlocksWorld:
    def __init__(self, stacks):
        self.stacks = stacks
        
    def get_top_blocks(self):
        return [stack[-1] if stack else None for stack in self.stacks]
    
    def is_valid_move(self, from_stack, to_stack):
        if from_stack < 0 or from_stack >= len(self.stacks) or not self.stacks[from_stack]:
            return False
        if to_stack < 0 or to_stack >= len(self.stacks):
            return False
        return True
    
    def make_move(self, from_stack, to_stack):
        if not self.is_valid_move(from_stack, to_stack):
            return None
        new_stacks = copy.deepcopy(self.stacks)
        block = new_stacks[from_stack].pop()
        new_stacks[to_stack].append(block)
        return BlocksWorld(new_stacks)
    
    def __str__(self):
        return str(self.stacks)
    
    def __eq__(self, other):
        return str(self) == str(other)
    
    def __hash__(self):
        return hash(str(self))

def find_solution(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = set()
    
    while queue:
        current_state, moves = queue.popleft()
        if current_state == goal_state:
            return moves
            
        if str(current_state) in visited:
            continue
            
        visited.add(str(current_state))
        
        for from_stack in range(len(current_state.stacks)):
            for to_stack in range(len(current_state.stacks)):
                if from_stack != to_stack and current_state.is_valid_move(from_stack, to_stack):
                    new_state = current_state.make_move(from_stack, to_stack)
                    if new_state:
                        block = current_state.stacks[from_stack][-1]
                        move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                        queue.append((new_state, moves + [move]))

# Initial and goal states
initial = BlocksWorld([['D'], ['B', 'E'], ['F', 'C', 'G', 'A']])
goal = BlocksWorld([['B', 'F', 'G'], ['A', 'C', 'D'], ['E']])

# Find solution
solution = find_solution(initial, goal)
if solution:
    print("<<<")
    print("\n".join(solution))
    print(">>>")
else:
    print("No solution found")