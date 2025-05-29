def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2, disk_moves):
    if n == 0:
        return
    if n == 1:
        disk_moves.append((1, source, target))
        return
    
    # Find the optimal k
    min_moves = float('inf')
    best_k = 0
    for k in range(1, n):
        moves = 2 * (3**k - 1) + (3**(n-k) - 1)
        if moves < min_moves:
            min_moves = moves
            best_k = k
    
    # Perform the moves
    hanoi_4_pegs(best_k, source, auxiliary1, auxiliary2, target, disk_moves)
    hanoi_3_pegs(n-best_k, source, target, auxiliary2, disk_moves)
    hanoi_4_pegs(best_k, auxiliary1, target, source, auxiliary2, disk_moves)

def hanoi_3_pegs(n, source, target, auxiliary, disk_moves):
    if n == 0:
        return
    hanoi_3_pegs(n-1, source, auxiliary, target, disk_moves)
    disk_moves.append((n, source, target))
    hanoi_3_pegs(n-1, auxiliary, target, source, disk_moves)

def format_moves(moves):
    return [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]

# Solve the problem
disk_moves = []
hanoi_4_pegs(6, 3, 4, 1, 2, disk_moves)
formatted_moves = format_moves(disk_moves)
print("\n".join(formatted_moves))