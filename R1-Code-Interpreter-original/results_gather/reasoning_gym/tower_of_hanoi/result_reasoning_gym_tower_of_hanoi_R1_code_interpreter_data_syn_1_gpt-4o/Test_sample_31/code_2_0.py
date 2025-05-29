def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2, start_disk=1):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk {start_disk} from Peg {source} to Peg {target}"]
    
    moves = []
    k = n // 2  # A common heuristic for choosing k
    # Move k disks to auxiliary1
    moves += hanoi_4_pegs(k, source, auxiliary1, target, auxiliary2, start_disk)
    # Move remaining n-k disks to target using auxiliary2
    moves += hanoi_4_pegs(n-k, source, target, auxiliary2, auxiliary1, start_disk + k)
    # Move k disks from auxiliary1 to target
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2, start_disk)
    
    return moves

# Solve the problem for 6 disks, moving from Peg 3 to Peg 2
moves = hanoi_4_pegs(6, 3, 2, 1, 4)
print("\n".join(moves))