from collections import deque
import copy

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
    # Get the block from top of source stack
    block = stacks[from_stack-1].pop()
    # Add it to destination stack
    stacks[to_stack-1].append(block)
    return block

def get_state_string(stacks):
    return "|".join(",".join(stack) for stack in stacks)

def find_solution():
    # Initial state
    initial_stacks = [
        [],                    # stack1
        ['G','A','E','F','C','B'],  # stack2
        ['D']                  # stack3
    ]
    
    # Goal state
    goal_stacks = [
        ['A','D','E'],        # stack1
        ['G'],                # stack2
        ['B','C','F']         # stack3
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
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack:
                    new_stacks = copy.deepcopy(current_stacks)
                    if is_valid_move(new_stacks, from_stack, to_stack):
                        block = make_move(new_stacks, from_stack, to_stack)
                        new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                        queue.append((new_stacks, new_moves))

solution = find_solution()
print("<<<" + "\n".join(solution) + ">>>")