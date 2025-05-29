from collections import deque
import copy

def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack has blocks
    if not stacks[from_stack]:
        return False
    return True

def apply_move(stacks, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
    new_stacks = copy.deepcopy(stacks)
    block = new_stacks[from_stack].pop()
    new_stacks[to_stack].append(block)
    return new_stacks

def state_to_string(stacks):
    return ';'.join(','.join(stack) for stack in stacks)

def solve_blocksworld():
    # Initial state
    stacks = [['F', 'A', 'G'], ['E', 'B'], ['C', 'D']]
    moves = []
    
    # Goal state
    goal = [['A', 'B', 'D'], ['G'], ['C', 'E', 'F']]
    
    # First, move G to stack2
    if stacks[0][-1] == 'G':
        moves.append("Move G from 1 to 2")
        stacks = apply_move(stacks, 0, 1)
    
    # Move A from stack1 to stack3
    if stacks[0][-1] == 'A':
        moves.append("Move A from 1 to 3")
        stacks = apply_move(stacks, 0, 2)
    
    # Move F from stack1 to stack3
    if stacks[0][-1] == 'F':
        moves.append("Move F from 1 to 3")
        stacks = apply_move(stacks, 0, 2)
    
    # Move E from stack2 to stack3
    if stacks[1] and stacks[1][-1] == 'E':
        moves.append("Move E from 2 to 3")
        stacks = apply_move(stacks, 1, 2)
    
    # Move B from stack2 to stack1
    if stacks[1] and stacks[1][-1] == 'B':
        moves.append("Move B from 2 to 1")
        stacks = apply_move(stacks, 1, 0)
    
    # Move D from stack3 to stack1
    if stacks[2][-1] == 'D':
        moves.append("Move D from 3 to 1")
        stacks = apply_move(stacks, 2, 0)

    print('\n'.join(moves))

solve_blocksworld()