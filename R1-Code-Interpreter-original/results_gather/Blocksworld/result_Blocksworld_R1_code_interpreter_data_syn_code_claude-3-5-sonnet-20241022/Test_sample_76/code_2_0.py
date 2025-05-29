from collections import deque
import copy

def is_valid_move(stacks, from_stack, to_stack):
    if not (0 <= from_stack-1 < len(stacks)) or not (0 <= to_stack-1 < len(stacks)):
        return False
    return len(stacks[from_stack-1]) > 0

def make_move(stacks, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
    new_stacks = copy.deepcopy(stacks)
    block = new_stacks[from_stack-1].pop()
    new_stacks[to_stack-1].append(block)
    return new_stacks

def evaluate_state(current, goal):
    score = 0
    # Check blocks in correct final position
    for i in range(len(current)):
        curr_stack = current[i]
        goal_stack = goal[i]
        
        # Penalty for blocks not in their final position
        correct_pos = 0
        for j in range(min(len(curr_stack), len(goal_stack))):
            if curr_stack[j] == goal_stack[j]:
                correct_pos += 1
            else:
                break
        score -= (len(curr_stack) - correct_pos) * 2
        
        # Extra penalty for blocks that need to be moved
        score -= abs(len(curr_stack) - len(goal_stack))
    
    return score

def find_solution():
    initial = [['A', 'F', 'D', 'E', 'H'], ['B', 'G'], ['C']]
    goal = [['G'], ['A', 'B', 'D', 'E'], ['C', 'F', 'H']]
    
    queue = [(initial, [])]
    visited = {str(initial)}
    best_score = float('-inf')
    best_solution = None
    
    while queue and (best_solution is None or len(queue[0][1]) < len(best_solution)):
        current_state, moves = queue.pop(0)
        current_score = evaluate_state(current_state, goal)
        
        if current_state == goal:
            if best_solution is None or len(moves) < len(best_solution):
                best_solution = moves
            continue
            
        if current_score > best_score:
            best_score = current_score
        
        # Generate and evaluate all possible moves
        possible_moves = []
        for from_stack in range(1, 4):
            if not current_state[from_stack-1]:
                continue
            moving_block = current_state[from_stack-1][-1]
            
            for to_stack in range(1, 4):
                if from_stack != to_stack:
                    new_state = make_move(current_state, from_stack, to_stack)
                    if new_state:
                        state_str = str(new_state)
                        if state_str not in visited:
                            new_moves = moves + [f"Move {moving_block} from {from_stack} to {to_stack}"]
                            score = evaluate_state(new_state, goal) - len(new_moves)
                            possible_moves.append((score, new_state, new_moves))
                            visited.add(state_str)
        
        # Add moves to queue, prioritizing better scores
        possible_moves.sort(key=lambda x: x[0], reverse=True)
        queue.extend((state, move) for _, state, move in possible_moves)
        queue.sort(key=lambda x: len(x[1]))  # Prioritize shorter solutions

    return best_solution

solution = find_solution()
if solution:
    print('\n'.join(solution))
else:
    print("No solution found")