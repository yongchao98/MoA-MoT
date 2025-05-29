def tower_of_hanoi(n, source, target, auxiliary, moves):
    if n == 1:
        moves.append(f"Move disk 1 from Peg {source} to Peg {target}")
        return
    tower_of_hanoi(n-1, source, auxiliary, target, moves)
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    tower_of_hanoi(n-1, auxiliary, target, source, moves)

# List to store the sequence of moves
moves = []

# Solve the Tower of Hanoi problem for 7 disks, moving from Peg 2 to Peg 3
tower_of_hanoi(7, 2, 3, 1, moves)

# Print all moves at once
for move in moves:
    print(move)