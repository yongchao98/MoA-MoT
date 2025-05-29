def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    # Find the optimal k
    min_moves = float('inf')
    best_k = 0
    for k in range(1, n):
        moves = 2 * len(hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)) + len(hanoi_4_pegs(n-k, source, target, auxiliary1, auxiliary2))
        if moves < min_moves:
            min_moves = moves
            best_k = k
    
    # Move the top k disks to auxiliary1
    moves = hanoi_4_pegs(best_k, source, auxiliary1, auxiliary2, target)
    # Move the remaining n-k disks to target
    moves += hanoi_4_pegs(n-best_k, source, target, auxiliary1, auxiliary2)
    # Move the k disks from auxiliary1 to target
    moves += hanoi_4_pegs(best_k, auxiliary1, target, source, auxiliary2)
    
    return moves

# Generate the sequence of moves
moves = hanoi_4_pegs(6, 3, 4, 1, 2)

# Format the output
formatted_moves = [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]
print("\n".join(formatted_moves))