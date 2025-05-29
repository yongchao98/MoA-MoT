def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return block

def print_solution(moves):
    result = "<<<" + "\n".join(moves) + ">>>"
    print(result)

# Initial state
stacks = [
    ['B', 'H', 'G', 'I', 'D'],  # stack 1
    ['A', 'F'],                 # stack 2
    ['E'],                      # stack 3
    ['J', 'C']                 # stack 4
]

# Goal state
goal = [
    ['H'],                      # stack 1
    ['A', 'B', 'D', 'G'],      # stack 2
    ['C', 'E', 'F', 'I', 'J']  # stack 3
]

# Store moves
moves = []

# Step by step solution
# Move D from stack 1 to stack 4
moves.append(f"Move D from 1 to 4")
make_move(stacks, 1, 4)

# Move I from stack 1 to stack 3
moves.append(f"Move I from 1 to 3")
make_move(stacks, 1, 3)

# Move G from stack 1 to stack 4
moves.append(f"Move G from 1 to 4")
make_move(stacks, 1, 4)

# Move H from stack 1 to stack 3
moves.append(f"Move H from 1 to 3")
make_move(stacks, 1, 3)

# Move B from stack 1 to stack 2
moves.append(f"Move B from 1 to 2")
make_move(stacks, 1, 2)

# Move F from stack 2 to stack 3
moves.append(f"Move F from 2 to 3")
make_move(stacks, 2, 3)

# Move G from stack 4 to stack 2
moves.append(f"Move G from 4 to 2")
make_move(stacks, 4, 2)

# Move D from stack 4 to stack 2
moves.append(f"Move D from 4 to 2")
make_move(stacks, 4, 2)

# Move C from stack 4 to stack 3
moves.append(f"Move C from 4 to 3")
make_move(stacks, 4, 3)

# Move J from stack 4 to stack 3
moves.append(f"Move J from 4 to 3")
make_move(stacks, 4, 3)

# Move H from stack 3 to stack 1
moves.append(f"Move H from 3 to 1")
make_move(stacks, 3, 1)

print_solution(moves)