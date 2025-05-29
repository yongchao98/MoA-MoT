def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2, moves):
    if n == 0:
        return
    if n == 1:
        moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
        return

    # Find the optimal k
    k = n // 2

    # Move k disks to auxiliary1
    hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target, moves)

    # Move n-k disks to target
    hanoi_4_pegs(n - k, source, target, auxiliary2, auxiliary1, moves)

    # Move k disks from auxiliary1 to target
    hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2, moves)

def move_disks(n, source, target, auxiliary1, auxiliary2):
    moves = []
    hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2, moves)
    return moves

# Generate the sequence of moves for 6 disks
moves = move_disks(6, 4, 1, 2, 3)
print("\n".join(moves))