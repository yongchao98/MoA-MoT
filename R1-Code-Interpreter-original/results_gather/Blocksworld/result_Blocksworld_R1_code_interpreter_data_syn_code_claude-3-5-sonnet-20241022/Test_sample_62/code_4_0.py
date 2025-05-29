def find_solution(initial, goal):
    moves = []
    
    # Helper function to get the top block of a stack
    def get_top(stack):
        return stack[-1] if stack else None
    
    # Helper function to find which stack contains a block
    def find_block(stacks, block):
        for i, stack in enumerate(stacks):
            if block in stack and stack[-1] == block:
                return i + 1
        return None
    
    current = [stack[:] for stack in initial]
    
    # First, move blocks to clear the way for C and G in stack2
    if 'H' in current[2]:  # Move H if it's on top
        moves.append(f"Move H from 3 to 1")
        current[0].append(current[2].pop())
    
    # Move C to stack2
    if 'C' in current[2]:
        moves.append(f"Move C from 3 to 2")
        current[1].append(current[2].pop())
    
    # Clear stack2 except for C
    while len(current[1]) > 1:
        block = current[1].pop()
        if block != 'C':
            moves.append(f"Move {block} from 2 to 3")
            current[2].append(block)
    
    # Move G to stack2
    if 'G' in current[1]:
        top = current[1].pop()
        while top != 'G':
            moves.append(f"Move {top} from 2 to 3")
            current[2].append(top)
            if current[1]:
                top = current[1].pop()
    moves.append(f"Move G from 2 to 2")
    current[1].append('G')
    
    # Move remaining blocks to stack3 in correct order
    blocks_for_stack3 = ['A', 'D', 'E', 'F', 'H']
    for block in blocks_for_stack3:
        if block in current[0]:
            moves.append(f"Move {block} from 1 to 3")
            current[2].append(current[0].pop())
        elif block in current[1] and block != 'G':
            moves.append(f"Move {block} from 2 to 3")
            current[2].append(current[1].pop())
    
    print("<<<")
    print("\n".join(moves))
    print(">>>")

# Initial state
initial = [['B'], ['G', 'A', 'F', 'E'], ['D', 'C', 'H']]
# Goal state
goal = [['B'], ['C', 'G'], ['A', 'D', 'E', 'F', 'H']]

find_solution(initial, goal)