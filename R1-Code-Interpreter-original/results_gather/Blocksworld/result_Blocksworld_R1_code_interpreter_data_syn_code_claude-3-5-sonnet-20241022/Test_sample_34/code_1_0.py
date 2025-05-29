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
    new_stacks = copy.deepcopy(stacks)
    block = new_stacks[from_stack-1].pop()
    new_stacks[to_stack-1].append(block)
    return new_stacks

def get_state_string(stacks):
    return ';'.join([''.join(stack) for stack in stacks])

def heuristic(current, goal):
    # Count blocks in wrong positions
    score = 0
    for curr_stack, goal_stack in zip(current, goal):
        curr_len = len(curr_stack)
        goal_len = len(goal_stack)
        for i in range(min(curr_len, goal_len)):
            if curr_stack[i] != goal_stack[i]:
                score += 1
    return score

def find_solution():
    initial_stacks = [['G', 'E', 'B'], ['A'], ['C', 'D', 'H', 'F']]
    goal_stacks = [['D', 'F', 'G'], ['A', 'C', 'E', 'H'], ['B']]
    
    # Priority queue with (heuristic_score + moves_made, current_stacks, moves)
    queue = [(heuristic(initial_stacks, goal_stacks), initial_stacks, [])]
    visited = {get_state_string(initial_stacks)}
    
    while queue:
        _, current_stacks, moves = queue.pop(0)
        
        if current_stacks == goal_stacks:
            return moves
        
        next_moves = []
        for from_stack in range(1, 4):
            for to_stack in range(1, 4):
                if from_stack != to_stack and current_stacks[from_stack-1]:
                    new_stacks = make_move(current_stacks, from_stack, to_stack)
                    if new_stacks:
                        state_string = get_state_string(new_stacks)
                        if state_string not in visited:
                            visited.add(state_string)
                            block = current_stacks[from_stack-1][-1]
                            new_moves = moves + [f"Move {block} from {from_stack} to {to_stack}"]
                            score = heuristic(new_stacks, goal_stacks) + len(new_moves)
                            next_moves.append((score, new_stacks, new_moves))
        
        # Add sorted next moves to queue
        queue.extend(sorted(next_moves))
    
    return None

solution = find_solution()
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found")