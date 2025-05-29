from collections import deque
import copy

def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack has blocks
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    new_stacks = copy.deepcopy(stacks)
    block = new_stacks[from_stack-1].pop()
    new_stacks[to_stack-1].append(block)
    return new_stacks

def is_block_in_goal_position(block, stack_idx, current_stacks, goal_stacks):
    # Check if block is in its final position in the goal state
    if stack_idx >= len(goal_stacks):
        return False
    current_stack = current_stacks[stack_idx]
    goal_stack = goal_stacks[stack_idx]
    if not current_stack:
        return False
    pos = len(current_stack) - 1
    return pos < len(goal_stack) and goal_stack[pos] == block

def find_solution(initial_stacks, goal_stacks):
    queue = deque([(initial_stacks, [])])
    visited = set()
    
    while queue:
        current_stacks, moves = queue.popleft()
        state_str = str(current_stacks)
        
        if current_stacks == goal_stacks:
            return moves
            
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        # Try all possible moves
        for from_stack in range(1, 4):
            if not current_stacks[from_stack-1]:  # Skip empty stacks
                continue
                
            block = current_stacks[from_stack-1][-1]
            
            for to_stack in range(1, 4):
                if from_stack == to_stack:
                    continue
                    
                if is_valid_move(current_stacks, from_stack, to_stack):
                    new_stacks = make_move(current_stacks, from_stack, to_stack)
                    
                    # Prioritize moves that put blocks in their goal position
                    move = f"Move {block} from {from_stack} to {to_stack}"
                    new_moves = moves + [move]
                    
                    # If this move puts a block in its goal position, prioritize it
                    if is_block_in_goal_position(block, to_stack-1, new_stacks, goal_stacks):
                        queue.appendleft((new_stacks, new_moves))
                    else:
                        queue.append((new_stacks, new_moves))
    
    return None

# Initial state
initial = [['B'], ['G', 'A', 'F', 'E'], ['D', 'C', 'H']]
# Goal state
goal = [['B'], ['C', 'G'], ['A', 'D', 'E', 'F', 'H']]

solution = find_solution(initial, goal)
if solution:
    print("<<<")
    print("\n".join(solution))
    print(">>>")
else:
    print("No solution found")