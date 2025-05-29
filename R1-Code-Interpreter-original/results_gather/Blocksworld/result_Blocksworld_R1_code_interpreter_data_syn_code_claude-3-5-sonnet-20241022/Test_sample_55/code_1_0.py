def get_stack_differences(current, goal):
    differences = 0
    for curr_stack, goal_stack in zip(current, goal):
        curr_set = set(curr_stack)
        goal_set = set(goal_stack)
        differences += len(curr_set.symmetric_difference(goal_set))
    return differences

def find_solution():
    # Initial state
    initial = [
        ['J'],
        ['C', 'K', 'F', 'E', 'B', 'I'],
        ['G', 'A'],
        ['H', 'D']
    ]
    
    # Goal state
    goal = [
        ['F', 'G', 'H', 'I', 'J'],
        ['C', 'D', 'E', 'K'],
        ['A', 'B'],
        []
    ]
    
    moves = []
    current = initial
    
    # Direct moves to build goal stacks
    while current != goal:
        # First, handle stack1 (F G H I J)
        if 'F' in current[1] and current[1][-1] == 'F':
            moves.append("Move F from 2 to 1")
            current[1].pop()
            current[0].append('F')
        elif 'G' in current[2] and current[2][-1] == 'G':
            moves.append("Move G from 3 to 1")
            current[2].pop()
            current[0].append('G')
        elif 'H' in current[3] and current[3][-1] == 'H':
            moves.append("Move H from 4 to 1")
            current[3].pop()
            current[0].append('H')
        elif 'I' in current[1] and current[1][-1] == 'I':
            moves.append("Move I from 2 to 1")
            current[1].pop()
            current[0].append('I')
        
        # Handle stack2 (C D E K)
        elif 'K' in current[1] and current[1][-1] == 'K':
            moves.append("Move K from 2 to 2")
            current[1].pop()
            current[1].append('K')
        elif 'D' in current[3] and current[3][-1] == 'D':
            moves.append("Move D from 4 to 2")
            current[3].pop()
            current[1].append('D')
        elif 'E' in current[1] and current[1][-1] == 'E':
            moves.append("Move E from 2 to 2")
            current[1].pop()
            current[1].append('E')
        
        # Handle stack3 (A B)
        elif 'B' in current[1] and current[1][-1] == 'B':
            moves.append("Move B from 2 to 3")
            current[1].pop()
            current[2].append('B')
        elif 'A' in current[2] and current[2][-1] == 'A':
            moves.append("Move A from 3 to 3")
            current[2].pop()
            current[2].append('A')
        
        # Clear obstacles
        elif len(current[1]) > 0:
            block = current[1][-1]
            moves.append(f"Move {block} from 2 to 4")
            current[1].pop()
            current[3].append(block)
        elif len(current[2]) > 0:
            block = current[2][-1]
            moves.append(f"Move {block} from 3 to 4")
            current[2].pop()
            current[3].append(block)
        
        if len(moves) > 30:  # Safety limit
            break
    
    return moves

solution = find_solution()
print("<<<" + "\n".join(solution) + ">>>")