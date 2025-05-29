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

# Step by step solution to reach the goal state
# First clear stack 1 except H
moves.append("Move D from 1 to 4")
make_move(stacks, 1, 4)

moves.append("Move I from 1 to 3")
make_move(stacks, 1, 3)

moves.append("Move G from 1 to 2")
make_move(stacks, 1, 2)

moves.append("Move B from 1 to 2")
make_move(stacks, 1, 2)

# Move blocks to their final positions in stack 3
moves.append("Move C from 4 to 3")
make_move(stacks, 4, 3)

moves.append("Move D from 4 to 2")
make_move(stacks, 4, 2)

moves.append("Move J from 4 to 3")
make_move(stacks, 4, 3)

moves.append("Move F from 2 to 3")
make_move(stacks, 2, 3)

moves.append("Move I from 3 to 3")
make_move(stacks, 3, 3)

print_solution(moves)