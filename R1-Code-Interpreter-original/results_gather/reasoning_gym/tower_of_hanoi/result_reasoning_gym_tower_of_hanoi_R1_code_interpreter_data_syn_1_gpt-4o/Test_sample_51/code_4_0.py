def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2, moves):
    if n == 0:
        return
    if n == 1:
        moves.append(f"Move disk 1 from Peg {source} to Peg {target}")
        return

    # Use a dynamic approach to find the optimal k
    k = n // 2

    # Step 1: Move the top k disks to auxiliary1
    hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target, moves)

    # Step 2: Move the remaining n-k disks to target
    hanoi_3_pegs(n - k, source, target, auxiliary2, moves)

    # Step 3: Move the k disks from auxiliary1 to target
    hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2, moves)

def hanoi_3_pegs(n, source, target, auxiliary, moves):
    if n == 0:
        return
    if n == 1:
        moves.append(f"Move disk 1 from Peg {source} to Peg {target}")
        return

    hanoi_3_pegs(n - 1, source, auxiliary, target, moves)
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    hanoi_3_pegs(n - 1, auxiliary, target, source, moves)

# Initialize the list to store moves
moves = []

# Solve the problem for 6 disks and 4 pegs
hanoi_4_pegs(6, 1, 4, 2, 3, moves)

# Print the sequence of moves
for move in moves:
    print(move)