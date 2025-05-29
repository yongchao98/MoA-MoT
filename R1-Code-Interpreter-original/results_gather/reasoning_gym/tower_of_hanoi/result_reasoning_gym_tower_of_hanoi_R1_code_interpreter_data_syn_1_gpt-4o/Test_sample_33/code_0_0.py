def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    # Find the optimal k
    min_moves = float('inf')
    best_k = 0
    for k in range(1, n):
        moves = 2 * hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target) + hanoi_4_pegs(n-k, source, target, auxiliary1, auxiliary2)
        if len(moves) < min_moves:
            min_moves = len(moves)
            best_k = k
    
    # Perform the moves
    moves = []
    moves += hanoi_4_pegs(best_k, source, auxiliary1, auxiliary2, target)
    moves += hanoi_4_pegs(n-best_k, source, target, auxiliary1, auxiliary2)
    moves += hanoi_4_pegs(best_k, auxiliary1, target, source, auxiliary2)
    
    return moves

def format_moves(moves):
    return [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]

# Solve the problem
moves = hanoi_4_pegs(6, 3, 4, 1, 2)
formatted_moves = format_moves(moves)
print("\n".join(formatted_moves))