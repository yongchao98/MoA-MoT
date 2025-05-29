def get_state_string(stacks):
    return ';'.join(','.join(s) for s in stacks)

def is_valid_move(stacks, from_idx, to_idx):
    return bool(stacks[from_idx])  # Source stack must have blocks

def make_move(stacks, from_idx, to_idx):
    block = stacks[from_idx].pop()
    stacks[to_idx].append(block)
    return block

def solve_blocksworld():
    # Initial and goal states
    initial = [['C'], ['B', 'E', 'A', 'J', 'F'], ['G', 'I'], ['D', 'H']]
    goal = [['E'], ['A', 'B', 'C', 'D', 'H'], ['F', 'G', 'I', 'J']]
    
    # Direct solving based on goal state analysis
    moves = []
    current = [list(stack) for stack in initial]
    
    # Step 1: Move F and J to stack 3 to build F,G,I,J
    if 'F' in current[1]:
        moves.append(f"Move F from 2 to 3")
        make_move(current, 1, 2)
    if 'J' in current[1]:
        moves.append(f"Move J from 2 to 3")
        make_move(current, 1, 2)
    
    # Step 2: Move E to stack 1
    if 'E' in current[1]:
        moves.append(f"Move E from 2 to 1")
        make_move(current, 1, 0)
    
    # Step 3: Build A,B,C,D,H in stack 2
    if 'H' in current[3]:
        moves.append(f"Move H from 4 to 2")
        make_move(current, 3, 1)
    if 'D' in current[3]:
        moves.append(f"Move D from 4 to 2")
        make_move(current, 3, 1)
    if 'C' in current[0]:
        moves.append(f"Move C from 1 to 2")
        make_move(current, 0, 1)
    
    print('\n'.join(moves))

solve_blocksworld()