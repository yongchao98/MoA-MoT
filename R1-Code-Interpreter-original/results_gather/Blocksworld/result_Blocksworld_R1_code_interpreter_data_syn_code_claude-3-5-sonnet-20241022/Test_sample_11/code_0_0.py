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
        [],                  # stack1
        ['E'],              # stack2
        ['B', 'A', 'F', 'J'],  # stack3
        ['C', 'G', 'D', 'H', 'I']   # stack4
    ]
    
    # Goal state
    goal = [
        ['C', 'D'],         # stack1
        ['E', 'F', 'I'],    # stack2
        ['A', 'B', 'G', 'H', 'J']  # stack3
    ]
    
    moves = []
    
    # Step 1: Move I from stack4
    moves.append(('I', 4, 2))
    stacks[3].pop()
    
    # Step 2: Move H from stack4
    moves.append(('H', 4, 3))
    stacks[3].pop()
    
    # Step 3: Move D from stack4
    moves.append(('D', 4, 1))
    stacks[3].pop()
    
    # Step 4: Move G from stack4
    moves.append(('G', 4, 3))
    stacks[3].pop()
    
    # Step 5: Move C from stack4
    moves.append(('C', 4, 1))
    stacks[3].pop()
    
    # Step 6: Move J from stack3
    moves.append(('J', 3, 3))
    stacks[2].pop()
    
    # Step 7: Move F from stack3
    moves.append(('F', 3, 2))
    stacks[2].pop()
    
    print(print_moves(moves))

solve_blocksworld()