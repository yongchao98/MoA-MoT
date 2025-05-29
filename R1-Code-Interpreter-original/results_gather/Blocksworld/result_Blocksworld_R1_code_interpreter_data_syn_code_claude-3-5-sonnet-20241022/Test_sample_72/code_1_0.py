def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return f"Move {block} from {from_stack} to {to_stack}"

def find_solution():
    # Initial state
    current_stacks = [['F'], ['B'], ['D', 'A', 'E', 'C']]
    goal_stacks = [['F'], ['C', 'E'], ['A', 'B', 'D']]
    
    moves = []
    
    # Corrected solution
    moves.append(make_move(current_stacks, 3, 2))  # Move C from 3 to 2
    moves.append(make_move(current_stacks, 3, 1))  # Move E from 3 to 1
    moves.append(make_move(current_stacks, 3, 1))  # Move A from 3 to 1
    moves.append(make_move(current_stacks, 2, 3))  # Move B from 2 to 3
    moves.append(make_move(current_stacks, 1, 3))  # Move A from 1 to 3
    moves.append(make_move(current_stacks, 1, 2))  # Move E from 1 to 2
    
    print("<<<" + "\n".join(moves) + ">>>")

find_solution()