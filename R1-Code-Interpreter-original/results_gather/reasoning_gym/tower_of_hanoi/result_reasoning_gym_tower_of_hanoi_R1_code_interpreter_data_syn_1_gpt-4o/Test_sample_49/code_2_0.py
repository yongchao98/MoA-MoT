def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2, disk_offset=0):
    if n == 0:
        return []
    if n == 1:
        return [(disk_offset + 1, source, target)]
    
    # Find the optimal k
    min_moves = float('inf')
    best_k = 0
    for k in range(1, n):
        moves1 = hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target, disk_offset)
        moves2 = hanoi_4_pegs(n-k, source, target, auxiliary1, auxiliary2, disk_offset + k)
        moves3 = hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2, disk_offset)
        total_moves = len(moves1) + len(moves2) + len(moves3)
        if total_moves < min_moves:
            min_moves = total_moves
            best_k = k
    
    # Perform the moves
    moves = []
    moves += hanoi_4_pegs(best_k, source, auxiliary1, auxiliary2, target, disk_offset)
    moves += hanoi_4_pegs(n-best_k, source, target, auxiliary1, auxiliary2, disk_offset + best_k)
    moves += hanoi_4_pegs(best_k, auxiliary1, target, source, auxiliary2, disk_offset)
    
    return moves

def format_moves(moves):
    return [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]

# Solve the problem
moves = hanoi_4_pegs(7, 3, 4, 1, 2)
formatted_moves = format_moves(moves)
print("\n".join(formatted_moves))