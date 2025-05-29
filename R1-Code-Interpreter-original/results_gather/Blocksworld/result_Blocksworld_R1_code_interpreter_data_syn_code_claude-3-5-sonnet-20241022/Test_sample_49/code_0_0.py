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
    return ';'.join([','.join(stack) for stack in stacks])

def is_goal_state(stacks, goal):
    return all(stack == goal[i] for i, stack in enumerate(stacks))

def find_solution():
    # Initial state
    initial_stacks = [
        [],                  # stack1
        ['E','G','H','D','A','C'],  # stack2
        [],                  # stack3
        ['B','F']           # stack4
    ]
    
    # Goal state
    goal_stacks = [
        ['A','C','F','G','H'],  # stack1
        [],                      # stack2
        ['B','D','E'],          # stack3
    ]
    
    queue = deque([(initial_stacks, [])])
    visited = {get_state_string(initial_stacks)}
    
    while queue:
        current_stacks, moves = queue.popleft()
        
        if len(moves) > 20:  # Limit depth to avoid infinite loops
            continue
            
        # Check if current state matches goal state (for the first 3 stacks)
        if all(current_stacks[i] == goal_stacks[i] for i in range(3)):
            return moves
            
        # Try all possible moves
        for from_stack in range(len(current_stacks)):
            if not current_stacks[from_stack]:
                continue
                
            for to_stack in range(len(current_stacks)):
                if from_stack == to_stack:
                    continue
                    
                new_stacks = make_move(current_stacks, from_stack, to_stack)
                if not new_stacks:
                    continue
                    
                state_string = get_state_string(new_stacks)
                if state_string in visited:
                    continue
                    
                block = current_stacks[from_stack][-1]
                move = f"Move {block} from {from_stack + 1} to {to_stack + 1}"
                new_moves = moves + [move]
                
                queue.append((new_stacks, new_moves))
                visited.add(state_string)
    
    return None

solution = find_solution()
if solution:
    print('\n'.join(solution))
else:
    print("No solution found")