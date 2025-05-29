from collections import deque
import copy

def get_block_goal_stack(block, goal_stacks):
    # Find which stack a block should end up in
    for i, stack in enumerate(goal_stacks):
        if block in stack:
            return i
    return None

def is_valid_move(stacks, from_stack, to_stack):
    if not stacks[from_stack-1]:
        return False
    return True

def evaluate_move(block, to_stack, current_stacks, goal_stacks):
    # Higher score is better
    score = 0
    goal_stack_idx = get_block_goal_stack(block, goal_stacks)
    
    # If moving to the correct stack
    if goal_stack_idx == to_stack - 1:
        goal_stack = goal_stacks[goal_stack_idx]
        current_pos = len(current_stacks[to_stack-1])
        if current_pos < len(goal_stack) and goal_stack[current_pos] == block:
            score += 20  # Highly prioritize correct positioning
    
    # Penalize moves that disturb correctly placed blocks
    if current_stacks[to_stack-1]:
        if to_stack-1 == get_block_goal_stack(current_stacks[to_stack-1][-1], goal_stacks):
            score -= 10
            
    return score

def find_solution(initial_stacks, goal_stacks):
    queue = []  # Priority queue
    visited = set()
    
    # Initialize with starting state
    queue.append((0, initial_stacks, []))
    
    while queue:
        _, current_stacks, moves = queue.pop(0)
        state_str = str(current_stacks)
        
        if current_stacks == goal_stacks:
            return moves
            
        if state_str in visited:
            continue
            
        visited.add(state_str)
        
        possible_moves = []
        for from_stack in range(1, 4):
            if not current_stacks[from_stack-1]:
                continue
                
            block = current_stacks[from_stack-1][-1]
            
            for to_stack in range(1, 4):
                if from_stack == to_stack:
                    continue
                    
                if is_valid_move(current_stacks, from_stack, to_stack):
                    new_stacks = copy.deepcopy(current_stacks)
                    new_stacks[from_stack-1].pop()
                    new_stacks[to_stack-1].append(block)
                    
                    score = evaluate_move(block, to_stack, current_stacks, goal_stacks)
                    move = f"Move {block} from {from_stack} to {to_stack}"
                    
                    # Priority is negative score (so higher scores are processed first)
                    priority = -score + len(moves)
                    possible_moves.append((priority, new_stacks, moves + [move]))
        
        # Add moves sorted by priority
        queue.extend(sorted(possible_moves))
    
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