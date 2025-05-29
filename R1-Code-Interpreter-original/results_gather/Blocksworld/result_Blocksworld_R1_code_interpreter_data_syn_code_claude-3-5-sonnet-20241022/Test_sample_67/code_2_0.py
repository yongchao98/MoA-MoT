def print_stacks(stacks):
    for i, stack in enumerate(stacks):
        print(f"Stack {i+1}: {' '.join(stack)}")

def is_valid_move(stacks, from_stack, to_stack):
    return len(stacks[from_stack-1]) > 0

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

moves = []

# Clear top of stack 1
moves.append("Move D from 1 to 4")
make_move(stacks, 1, 4)

moves.append("Move I from 1 to 3")
make_move(stacks, 1, 3)

moves.append("Move G from 1 to 4")
make_move(stacks, 1, 4)

# Move B to its final position
moves.append("Move H from 1 to 3")
make_move(stacks, 1, 3)

moves.append("Move B from 1 to 2")
make_move(stacks, 1, 2)

# Move F to stack 3
moves.append("Move F from 2 to 3")
make_move(stacks, 2, 3)

# Arrange stack 2
moves.append("Move G from 4 to 2")
make_move(stacks, 4, 2)

moves.append("Move D from 4 to 2")
make_move(stacks, 4, 2)

# Complete stack 3
moves.append("Move C from 4 to 3")
make_move(stacks, 4, 3)

moves.append("Move J from 4 to 3")
make_move(stacks, 4, 3)

# Finally place H in stack 1
moves.append("Move H from 3 to 1")
make_move(stacks, 3, 1)

print_solution(moves)