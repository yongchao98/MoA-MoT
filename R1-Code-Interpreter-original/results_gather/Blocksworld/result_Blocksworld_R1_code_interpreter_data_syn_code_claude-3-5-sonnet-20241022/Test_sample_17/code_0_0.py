class State:
    def __init__(self, stacks):
        self.stacks = stacks
        
    def is_valid_move(self, from_stack, to_stack):
        # Check if source stack is not empty
        if not self.stacks[from_stack]:
            return False
        return True
        
    def move(self, from_stack, to_stack):
        if self.is_valid_move(from_stack, to_stack):
            block = self.stacks[from_stack][-1]
            new_stacks = [stack[:] for stack in self.stacks]
            new_stacks[from_stack].pop()
            new_stacks[to_stack].append(block)
            return State(new_stacks), f"Move {block} from {from_stack+1} to {to_stack+1}"
        return None, None

    def __eq__(self, other):
        return all(self.stacks[i] == other.stacks[i] for i in range(len(self.stacks)))
        
    def __hash__(self):
        return hash(tuple(tuple(stack) for stack in self.stacks))

def find_solution(initial_state, goal_state):
    from collections import deque
    
    queue = deque([(initial_state, [])])
    visited = set()
    
    while queue:
        current_state, moves = queue.popleft()
        if current_state == goal_state:
            return moves
            
        state_hash = hash(current_state)
        if state_hash in visited:
            continue
        visited.add(state_hash)
        
        for from_stack in range(len(current_state.stacks)):
            for to_stack in range(len(current_state.stacks)):
                if from_stack != to_stack:
                    next_state, move = current_state.move(from_stack, to_stack)
                    if next_state:
                        queue.append((next_state, moves + [move]))
    return None

# Initial and goal states
initial = State([
    ['C', 'F', 'E'],  # stack 1
    ['D'],            # stack 2
    [],               # stack 3
    ['A', 'B']        # stack 4
])

goal = State([
    ['C', 'E', 'F'],  # stack 1
    ['B', 'D'],       # stack 2
    ['A'],            # stack 3
    []                # stack 4
])

# Find and print solution
solution = find_solution(initial, goal)
if solution:
    print("<<<")
    print("\n".join(solution))
    print(">>>")
else:
    print("No solution found")