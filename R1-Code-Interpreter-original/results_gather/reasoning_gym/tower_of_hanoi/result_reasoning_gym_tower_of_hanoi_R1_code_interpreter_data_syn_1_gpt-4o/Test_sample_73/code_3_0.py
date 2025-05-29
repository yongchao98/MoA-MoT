def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2, disk_num):
    if n == 0:
        return []
    if n == 1:
        return [(disk_num, source, target)]
    
    # Find the optimal k
    min_moves = float('inf')
    best_k = 0
    for k in range(1, n):
        moves1 = hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target, disk_num)
        moves2 = hanoi_3_pegs(n - k, source, target, auxiliary2, disk_num + k)
        moves3 = hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2, disk_num)
        total_moves = len(moves1) + len(moves2) + len(moves3)
        if total_moves < min_moves:
            min_moves = total_moves
            best_k = k
    
    # Perform the moves
    moves = []
    moves += hanoi_4_pegs(best_k, source, auxiliary1, auxiliary2, target, disk_num)
    moves += hanoi_3_pegs(n - best_k, source, target, auxiliary2, disk_num + best_k)
    moves += hanoi_4_pegs(best_k, auxiliary1, target, source, auxiliary2, disk_num)
    
    return moves

def hanoi_3_pegs(n, source, target, auxiliary, disk_num):
    if n == 0:
        return []
    if n == 1:
        return [(disk_num, source, target)]
    moves = []
    moves += hanoi_3_pegs(n - 1, source, auxiliary, target, disk_num)
    moves.append((disk_num + n - 1, source, target))
    moves += hanoi_3_pegs(n - 1, auxiliary, target, source, disk_num)
    return moves

def format_moves(moves):
    return [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]

# Solve the problem
moves = hanoi_4_pegs(7, 1, 4, 2, 3, 1)
formatted_moves = format_moves(moves)
print("\n".join(formatted_moves))