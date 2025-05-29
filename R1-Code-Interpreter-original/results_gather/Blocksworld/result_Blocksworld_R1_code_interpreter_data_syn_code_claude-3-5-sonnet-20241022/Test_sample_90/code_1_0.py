def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(stacks, from_stack, to_stack):
    if not stacks[from_stack]:  # Source stack is empty
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
    block = stacks[from_stack].pop()
    stacks[to_stack].append(block)
    return block

def solve_blocksworld():
    # Initial and goal states
    initial = [['D', 'E', 'G'], ['C', 'B', 'A'], ['F', 'H']]
    goal = [['B', 'E', 'G'], ['C', 'D', 'H'], ['A', 'F']]
    
    moves = []
    current = [list(stack) for stack in initial]
    
    # First, move A to stack 3
    if 'A' in current[1]:
        moves.append("Move A from 2 to 3")
        make_move(current, 1, 2)
    
    # Move B to stack 1
    if 'B' in current[1]:
        moves.append("Move B from 2 to 1")
        make_move(current, 1, 0)
    
    # Move H to stack 2
    if 'H' in current[2]:
        moves.append("Move H from 3 to 2")
        make_move(current, 2, 1)
    
    # Move D to stack 2
    if 'D' in current[0]:
        moves.append("Move D from 1 to 2")
        make_move(current, 0, 1)

    print("<<<" + "\n".join(moves) + ">>>")

solve_blocksworld()