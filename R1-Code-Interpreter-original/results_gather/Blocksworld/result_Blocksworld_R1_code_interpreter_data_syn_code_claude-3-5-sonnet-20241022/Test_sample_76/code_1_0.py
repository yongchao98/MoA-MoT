from collections import deque
import copy

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty and stacks exist
    if not (0 <= from_stack-1 < len(stacks)) or not (0 <= to_stack-1 < len(stacks)):
        return False
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
    new_stacks = copy.deepcopy(stacks)
    block = new_stacks[from_stack-1].pop()
    new_stacks[to_stack-1].append(block)
    return new_stacks

def get_stack_differences(current, goal):
    differences = 0
    for curr_stack, goal_stack in zip(current, goal):
        differences += sum(1 for i in range(min(len(curr_stack), len(goal_stack))) if curr_stack[i] != goal_stack[i])
        differences += abs(len(curr_stack) - len(goal_stack))
    return differences

def find_solution():
    initial = [['A', 'F', 'D', 'E', 'H'], ['B', 'G'], ['C']]
    goal = [['G'], ['A', 'B', 'D', 'E'], ['C', 'F', 'H']]
    
    # Priority queue based on number of moves + heuristic
    queue = [(initial, [])]
    visited = {str(initial)}
    
    while queue:
        current_state, moves = queue.pop(0)
        
        if current_state == goal:
            return moves
            
        # Generate all possible moves
        possible_moves = []
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack and len(current_state[from_stack-1]) > 0:
                    new_state = make_move(current_state, from_stack, to_stack)
                    if new_state:
                        state_str = str(new_state)
                        if state_str not in visited:
                            block = current_state[from_stack-1][-1]
                            new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                            score = len(new_moves) + get_stack_differences(new_state, goal)
                            possible_moves.append((score, new_state, new_moves))
                            visited.add(state_str)
        
        # Sort moves by score (lower is better)
        possible_moves.sort(key=lambda x: x[0])
        queue.extend((state, move) for _, state, move in possible_moves)

    return None

solution = find_solution()
if solution:
    print('\n'.join(solution))
else:
    print("No solution found")