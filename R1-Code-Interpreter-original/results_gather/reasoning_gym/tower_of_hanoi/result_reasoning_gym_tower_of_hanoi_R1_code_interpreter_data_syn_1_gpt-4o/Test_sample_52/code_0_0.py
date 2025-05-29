def hanoi(n, source, target, auxiliary, moves):
    if n == 1:
        moves.append(f"Move disk 1 from Peg {source} to Peg {target}")
    else:
        hanoi(n-1, source, auxiliary, target, moves)
        moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
        hanoi(n-1, auxiliary, target, source, moves)

# Initialize parameters
n_disks = 5
source_peg = 2
target_peg = 1
auxiliary_peg = 3
moves = []

# Generate the sequence of moves
hanoi(n_disks, source_peg, target_peg, auxiliary_peg, moves)

# Print the sequence of moves
for move in moves:
    print(move)