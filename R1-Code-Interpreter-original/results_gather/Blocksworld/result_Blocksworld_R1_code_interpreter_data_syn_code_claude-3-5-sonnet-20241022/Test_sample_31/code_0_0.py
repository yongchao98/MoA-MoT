from collections import deque
import copy

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty and destination exists
    if not stacks[from_stack] or to_stack >= len(stacks):
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
    new_stacks = copy.deepcopy(stacks)
    block = new_stacks[from_stack].pop()
    new_stacks[to_stack].append(block)
    return new_stacks

def get_state_string(stacks):
    return ';'.join(','.join(stack) for stack in stacks)

def find_solution():
    # Initial state
    initial_stacks = [
        [],             # stack1
        ['B', 'E'],    # stack2
        ['C', 'F', 'D'],# stack3
        ['A']          # stack4
    ]
    
    # Goal state
    goal_stacks = [
        ['A', 'B', 'F'],# stack1
        ['C'],         # stack2
        ['D', 'E'],    # stack3
        []             # stack4
    ]
    
    # BFS search
    queue = deque([(initial_stacks, [])])
    visited = set()
    
    while queue:
        current_stacks, moves = queue.popleft()
        current_state = get_state_string(current_stacks)
        
        if current_state == get_state_string(goal_stacks):
            return moves
        
        if current_state in visited:
            continue
            
        visited.add(current_state)
        
        # Try all possible moves
        for from_stack in range(len(current_stacks)):
            for to_stack in range(len(current_stacks)):
                if from_stack != to_stack and current_stacks[from_stack]:
                    new_stacks = make_move(current_stacks, from_stack, to_stack)
                    if new_stacks:
                        block = current_stacks[from_stack][-1]
                        move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                        queue.append((new_stacks, moves + [move]))

solution = find_solution()
if solution:
    print('\n'.join(solution))