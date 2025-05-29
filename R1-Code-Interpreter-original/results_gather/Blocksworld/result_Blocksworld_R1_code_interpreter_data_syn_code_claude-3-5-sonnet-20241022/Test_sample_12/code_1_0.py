def print_moves(stacks, moves):
    print("<<<")
    for move in moves:
        print(move)
    print(">>>")

# Initial state
stacks = [
    [],  # stack 1
    ['K', 'D', 'I', 'G', 'B', 'J'],  # stack 2
    ['H', 'C', 'E', 'A', 'F']  # stack 3
]

# Record moves
moves = []

# First, clear the way to get F, E, A for stack 1
moves.append("Move J from 2 to 3")  # Move J to its final position in stack 3
moves.append("Move B from 2 to 1")  # Temporarily move B
moves.append("Move G from 2 to 2")  # Move G to its final position in stack 2
moves.append("Move I from 2 to 3")  # Move I to its final position in stack 3
moves.append("Move D from 2 to 3")  # Move D to its final position in stack 3
moves.append("Move K from 2 to 2")  # Move K to its final position in stack 2

# Now move F, A, E to stack 1
moves.append("Move F from 3 to 1")
moves.append("Move A from 3 to 1")
moves.append("Move E from 3 to 1")

# Finally, arrange remaining blocks in stack 3
moves.append("Move C from 3 to 3")  # C is already in position
moves.append("Move H from 3 to 3")  # H is already in position

print_moves(stacks, moves)