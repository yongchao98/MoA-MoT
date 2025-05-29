from collections import deque
import copy

def is_block_in_final_position(block, stack_idx, position, goal_stacks):
    # Check if this position for the block matches the goal state
    return (stack_idx < len(goal_stacks) and 
            position < len(goal_stacks[stack_idx]) and 
            goal_stacks[stack_idx][position] == block)

def evaluate_move(current_stacks, from_stack, to_stack, block, goal_stacks):
    # Higher score is better
    score = 0
    # Prefer moves that put blocks in their final position
    if any(goal_stack and goal_stack[-1] == block for goal_stack in goal_stacks):
        stack_idx = next(i for i, stack in enumerate(goal_stacks) if stack and stack[-1] == block)
        if stack_idx == to_stack - 1:
            score += 10
    return score

def find_solution(initial_stacks, goal_stacks):
    visited = set()
    # Priority queue with (priority, stacks, moves)
    queue = [(0, initial_stacks, [])]
    max_moves = 10  # Limit solution length
    
    while queue:
        priority, current_stacks, moves = queue.pop(0)
        
        if len(moves) > max_moves:
            continue
            
        current_state = str(current_stacks)
        if current_state in visited:
            continue
        visited.add(current_state)
        
        if current_stacks == goal_stacks:
            return moves
            
        possible_moves = []
        for from_stack in range(1, 4):
            if not current_stacks[from_stack-1]:  # Skip empty stacks
                continue
            block = current_stacks[from_stack-1][-1]
            
            for to_stack in range(1, 4):
                if from_stack != to_stack:
                    new_stacks = copy.deepcopy(current_stacks)
                    block = new_stacks[from_stack-1].pop()
                    new_stacks[to_stack-1].append(block)
                    
                    score = evaluate_move(current_stacks, from_stack, to_stack, block, goal_stacks)
                    move = f"Move {block} from {from_stack} to {to_stack}"
                    
                    # Priority is negative score (so lower scores are popped first)
                    # plus number of moves to prefer shorter solutions
                    priority = -score + len(moves)
                    possible_moves.append((priority, new_stacks, moves + [move]))
        
        # Add all possible moves, sorted by priority
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