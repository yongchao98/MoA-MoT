def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    # Find the optimal k
    min_moves = float('inf')
    best_k = 0
    for k in range(1, n):
        moves = 2 * len(hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)) + len(hanoi_3_pegs(n-k, source, target, auxiliary1))
        if moves < min_moves:
            min_moves = moves
            best_k = k
    
    # Move the top best_k disks to auxiliary1
    moves = hanoi_4_pegs(best_k, source, auxiliary1, auxiliary2, target)
    # Move the remaining n-best_k disks to target
    moves += hanoi_3_pegs(n-best_k, source, target, auxiliary2)
    # Move the best_k disks from auxiliary1 to target
    moves += hanoi_4_pegs(best_k, auxiliary1, target, source, auxiliary2)
    
    return moves

def hanoi_3_pegs(n, source, target, auxiliary):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    moves = hanoi_3_pegs(n-1, source, auxiliary, target)
    moves.append((n, source, target))
    moves += hanoi_3_pegs(n-1, auxiliary, target, source)
    
    return moves

def format_moves(moves):
    return [f"Move disk {disk} from Peg {src} to Peg {tgt}" for disk, src, tgt in moves]

# Solve the problem
moves = hanoi_4_pegs(7, 2, 3, 1, 4)
formatted_moves = format_moves(moves)
print("\n".join(formatted_moves))