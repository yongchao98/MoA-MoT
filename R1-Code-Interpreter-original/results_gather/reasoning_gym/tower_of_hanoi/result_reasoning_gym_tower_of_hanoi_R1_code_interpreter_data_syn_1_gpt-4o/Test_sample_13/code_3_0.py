def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(n, source, target)]
    
    # Find the optimal k
    min_moves = float('inf')
    best_k = 0
    for k in range(1, n):
        moves1 = hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
        moves2 = hanoi_4_pegs(n - k, source, target, auxiliary1, auxiliary2)
        moves3 = hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
        total_moves = len(moves1) + len(moves2) + len(moves3)
        if total_moves < min_moves:
            min_moves = total_moves
            best_k = k
    
    # Move the top k disks to auxiliary1
    moves = hanoi_4_pegs(best_k, source, auxiliary1, auxiliary2, target)
    # Move the remaining n-k disks to target
    moves += hanoi_4_pegs(n - best_k, source, target, auxiliary1, auxiliary2)
    # Move the k disks from auxiliary1 to target
    moves += hanoi_4_pegs(best_k, auxiliary1, target, source, auxiliary2)
    
    return moves

def format_moves(moves):
    return [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]

# Solve the problem for 5 disks
moves = hanoi_4_pegs(5, 4, 2, 1, 3)
formatted_moves = format_moves(moves)
print("\n".join(formatted_moves))