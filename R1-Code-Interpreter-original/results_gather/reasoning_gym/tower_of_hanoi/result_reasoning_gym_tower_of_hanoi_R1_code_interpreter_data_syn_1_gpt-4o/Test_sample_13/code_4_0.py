def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2, disks):
    if n == 0:
        return []
    if n == 1:
        return [(disks[0], source, target)]
    
    # Find the optimal k
    min_moves = float('inf')
    best_k = 0
    for k in range(1, n):
        moves1 = hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target, disks[:k])
        moves2 = hanoi_4_pegs(n - k, source, target, auxiliary1, auxiliary2, disks[k:])
        moves3 = hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2, disks[:k])
        total_moves = len(moves1) + len(moves2) + len(moves3)
        if total_moves < min_moves:
            min_moves = total_moves
            best_k = k
    
    # Move the top k disks to auxiliary1
    moves = hanoi_4_pegs(best_k, source, auxiliary1, auxiliary2, target, disks[:best_k])
    # Move the remaining n-k disks to target
    moves += hanoi_4_pegs(n - best_k, source, target, auxiliary1, auxiliary2, disks[best_k:])
    # Move the k disks from auxiliary1 to target
    moves += hanoi_4_pegs(best_k, auxiliary1, target, source, auxiliary2, disks[:best_k])
    
    return moves

def format_moves(moves):
    return [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]

# Solve the problem for 5 disks
disks = list(range(1, 6))  # Disk numbers 1 to 5
moves = hanoi_4_pegs(5, 4, 2, 1, 3, disks)
formatted_moves = format_moves(moves)
print("\n".join(formatted_moves))