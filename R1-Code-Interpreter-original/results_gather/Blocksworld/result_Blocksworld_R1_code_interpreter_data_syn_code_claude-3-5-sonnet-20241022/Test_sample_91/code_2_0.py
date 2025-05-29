def get_state_string(stacks):
    return ';'.join(','.join(stack) for stack in stacks)

def is_valid_state(state, goal):
    # Check if current state could potentially lead to goal
    blocks_count = {}
    for stack in state:
        for block in stack:
            blocks_count[block] = blocks_count.get(block, 0) + 1
    
    goal_blocks = {}
    for stack in goal:
        for block in stack:
            goal_blocks[block] = goal_blocks.get(block, 0) + 1
    
    return blocks_count == goal_blocks

def get_moves():
    # Initial and goal states
    initial = [['E'], ['B'], ['H', 'C', 'A', 'F'], ['D', 'G']]
    goal = [['D', 'E', 'G', 'H'], ['B', 'C', 'F'], ['A'], []]
    
    moves_list = []
    current = [stack[:] for stack in initial]
    
    # Direct solution based on problem analysis
    # Move F out of the way
    if 'F' in current[2]:
        moves_list.append("Move F from 3 to 4")
        current[2].remove('F')
        current[3].append('F')
    
    # Move A to its final position
    if 'A' in current[2]:
        moves_list.append("Move A from 3 to 3")
        current[2].remove('A')
        current[2].append('A')
    
    # Move C to stack 2
    if 'C' in current[2]:
        moves_list.append("Move C from 3 to 2")
        current[2].remove('C')
        current[1].append('C')
    
    # Move H to stack 1
    if 'H' in current[2]:
        moves_list.append("Move H from 3 to 1")
        current[2].remove('H')
        current[0].append('H')
    
    # Move G to stack 1
    if 'G' in current[3]:
        moves_list.append("Move G from 4 to 1")
        current[3].remove('G')
        current[0].append('G')
    
    # Move F to stack 2
    if 'F' in current[3]:
        moves_list.append("Move F from 4 to 2")
        current[3].remove('F')
        current[1].append('F')
    
    # Move D to stack 1
    if 'D' in current[3]:
        moves_list.append("Move D from 4 to 1")
        current[3].remove('D')
        current[0].append('D')
    
    # Move E to correct position in stack 1
    if 'E' in current[0]:
        moves_list.append("Move E from 1 to 1")
    
    return moves_list

# Run solution and print result
solution = get_moves()
print("<<<" + "\n".join(solution) + ">>>")