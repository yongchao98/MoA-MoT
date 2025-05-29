def move_disks(n, source, target, auxiliary1, auxiliary2, moves):
    if n == 0:
        return
    if n == 1:
        moves.append(f"Move disk 1 from Peg {source} to Peg {target}")
        return
    # Move top n-1 disks from source to auxiliary1
    move_disks(n-1, source, auxiliary1, auxiliary2, target, moves)
    # Move the nth disk from source to target
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    # Move the n-1 disks from auxiliary1 to target
    move_disks(n-1, auxiliary1, target, source, auxiliary2, moves)

# Initialize the list to store the moves
moves = []
# Solve the problem for 3 disks, moving from Peg 2 to Peg 4
move_disks(2, 2, 1, 3, 4, moves)  # Move top 2 disks to Peg 1
moves.append("Move disk 3 from Peg 2 to Peg 4")  # Move largest disk to Peg 4
move_disks(2, 1, 4, 2, 3, moves)  # Move 2 disks from Peg 1 to Peg 4

# Print the sequence of moves
for move in moves:
    print(move)