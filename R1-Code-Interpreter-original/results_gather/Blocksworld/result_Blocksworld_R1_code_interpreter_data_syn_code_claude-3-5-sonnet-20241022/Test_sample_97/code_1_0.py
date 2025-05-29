def get_heuristic(stacks, goal_stacks):
    # Simple heuristic: count blocks in wrong positions
    score = 0
    for i, stack in enumerate(stacks):
        goal_stack = goal_stacks[i]
        for j, block in enumerate(stack):
            if j < len(goal_stack) and block != goal_stack[j]:
                score += 1
            elif j >= len(goal_stack):
                score += 1
    return score

def solve_blocksworld():
    initial = [
        ['H', 'J', 'A', 'B'],
        ['F', 'I', 'K', 'C'],
        [],
        ['G', 'E', 'D']
    ]
    goal = [
        ['B', 'C', 'D', 'E', 'I', 'K'],
        ['A', 'H', 'J'],
        ['F', 'G'],
        []
    ]
    
    moves = []
    current = initial
    
    # Direct solving approach with specific strategy
    while current != goal:
        # Move D to position
        if 'D' in current[3]:
            moves.append("Move D from 4 to 1")
            current[1].append(current[3].pop())
            
        # Move C to position
        if 'C' in current[1]:
            moves.append("Move C from 2 to 1")
            current[0].append(current[1].pop())
            
        # Move B to position
        if 'B' in current[0] and current[0][-1] == 'B':
            moves.append("Move B from 1 to 3")
            current[2].append(current[0].pop())
            
        # Move A for H,J stack
        if 'A' in current[0] and current[0][-1] == 'A':
            moves.append("Move A from 1 to 3")
            current[2].append(current[0].pop())
            
        # Move J for H,J stack
        if 'J' in current[0] and current[0][-1] == 'J':
            moves.append("Move J from 1 to 2")
            current[1].append(current[0].pop())
            
        # Move H for H,J stack
        if 'H' in current[0] and current[0][-1] == 'H':
            moves.append("Move H from 1 to 2")
            current[1].append(current[0].pop())
            
        # Move K to position
        if 'K' in current[1] and current[1][-1] == 'K':
            moves.append("Move K from 2 to 1")
            current[0].append(current[1].pop())
            
        # Move I to position
        if 'I' in current[1] and current[1][-1] == 'I':
            moves.append("Move I from 2 to 1")
            current[0].append(current[1].pop())
            
        # Move F,G to position
        if 'F' in current[1] and current[1][-1] == 'F':
            moves.append("Move F from 2 to 3")
            current[2].append(current[1].pop())
            
        if 'G' in current[3] and current[3][-1] == 'G':
            moves.append("Move G from 4 to 3")
            current[2].append(current[3].pop())
            
        # Move E to position
        if 'E' in current[3] and current[3][-1] == 'E':
            moves.append("Move E from 4 to 1")
            current[0].append(current[3].pop())
            
        # Final adjustments
        if 'A' in current[2] and current[2][-1] == 'A':
            moves.append("Move A from 3 to 2")
            current[1].append(current[2].pop())
            
        if 'B' in current[2] and current[2][-1] == 'B':
            moves.append("Move B from 3 to 1")
            current[0].append(current[2].pop())
            
        break  # Prevent infinite loop
        
    print('\n'.join(moves))

solve_blocksworld()