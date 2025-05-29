def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def print_moves(moves):
    result = ""
    for move in moves:
        result += f"Move {move[0]} from {move[1]} to {move[2]}\n"
    return result.strip()

def get_top_block(stack):
    return stack[-1] if stack else None

def solve_blocksworld():
    # Initial state
    stacks = [
        ['D', 'E'],        # stack1
        ['A', 'I', 'G'],   # stack2
        ['C', 'B'],        # stack3
        ['K', 'F', 'J', 'H']  # stack4
    ]
    
    # Goal state
    goal = [
        ['A', 'C', 'D', 'G'],        # stack1
        ['E', 'F', 'H', 'I', 'J', 'K'],   # stack2
        ['B']         # stack3
    ]
    
    moves = []
    
    # Helper function to find a block's location
    def find_block(block):
        for i, stack in enumerate(stacks):
            if block in stack:
                return i + 1
        return None
    
    # Move H from stack4 to stack2 (it's on top)
    moves.append(('H', 4, 2))
    stacks[3].pop()
    stacks[1].append('H')
    
    # Move J from stack4 to stack2
    moves.append(('J', 4, 2))
    stacks[3].pop()
    stacks[1].append('J')
    
    # Move F from stack4 to stack2
    moves.append(('F', 4, 2))
    stacks[3].pop()
    stacks[1].append('F')
    
    # Move K from stack4 to stack2
    moves.append(('K', 4, 2))
    stacks[3].pop()
    stacks[1].append('K')
    
    # Move B from stack3 to stack3 (goal position)
    moves.append(('B', 3, 3))
    stacks[2].pop()
    stacks[2].append('B')
    
    # Move C from stack3 to stack1
    moves.append(('C', 3, 1))
    stacks[2].pop()
    stacks[1].append('C')
    
    # Move G from stack2 to stack1
    moves.append(('G', 2, 1))
    stacks[1].pop()
    stacks[0].append('G')
    
    # Move I from stack2 to stack2
    moves.append(('I', 2, 2))
    stacks[1].pop()
    stacks[1].append('I')
    
    # Move A from stack2 to stack1
    moves.append(('A', 2, 1))
    stacks[1].pop()
    stacks[0].append('A')
    
    # Move E from stack1 to stack2
    moves.append(('E', 1, 2))
    stacks[0].pop()
    stacks[1].append('E')
    
    # Move D from stack1 to stack1
    moves.append(('D', 1, 1))
    stacks[0].pop()
    stacks[0].append('D')
    
    print(print_moves(moves))

solve_blocksworld()