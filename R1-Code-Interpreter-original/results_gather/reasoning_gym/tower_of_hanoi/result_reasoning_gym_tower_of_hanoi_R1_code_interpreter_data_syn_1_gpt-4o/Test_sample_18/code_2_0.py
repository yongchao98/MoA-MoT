def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2, disks):
    if n == 0:
        return []
    if n == 1:
        return [(disks[0], source, target)]
    
    # Find the optimal k
    min_moves = float('inf')
    best_k = 0
    for k in range(1, n):
        moves = 2 * len(hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target, disks[:k])) + len(hanoi_4_pegs(n-k, source, target, auxiliary1, auxiliary2, disks[k:]))
        if moves < min_moves:
            min_moves = moves
            best_k = k
    
    # Move the top k disks to auxiliary1
    moves = hanoi_4_pegs(best_k, source, auxiliary1, auxiliary2, target, disks[:best_k])
    # Move the remaining n-k disks to target
    moves += hanoi_4_pegs(n-best_k, source, target, auxiliary1, auxiliary2, disks[best_k:])
    # Move the k disks from auxiliary1 to target
    moves += hanoi_4_pegs(best_k, auxiliary1, target, source, auxiliary2, disks[:best_k])
    
    return moves

# Generate the sequence of moves
disks = list(range(1, 7))  # Disks are numbered 1 to 6
moves = hanoi_4_pegs(6, 3, 4, 1, 2, disks)

# Format the output
formatted_moves = [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]
print("\n".join(formatted_moves))