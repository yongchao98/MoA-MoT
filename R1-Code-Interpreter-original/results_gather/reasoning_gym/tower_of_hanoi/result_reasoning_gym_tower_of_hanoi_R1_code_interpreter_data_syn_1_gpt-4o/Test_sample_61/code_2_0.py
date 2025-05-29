def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    moves = []
    k = 3  # For 6 disks, the optimal k is known to be 3
    # Move k disks to auxiliary1
    moves += hanoi_4_pegs(k, source, auxiliary1, target, auxiliary2)
    # Move n-k disks to target
    moves += hanoi_3_pegs(n-k, source, target, auxiliary2)
    # Move k disks from auxiliary1 to target
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    
    return moves

def hanoi_3_pegs(n, source, target, auxiliary):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    moves = []
    moves += hanoi_3_pegs(n-1, source, auxiliary, target)
    moves.append((n, source, target))
    moves += hanoi_3_pegs(n-1, auxiliary, target, source)
    
    return moves

def format_moves(moves):
    return [f"Move disk {disk} from Peg {src} to Peg {tgt}" for disk, src, tgt in moves]

# Solve the problem
moves = hanoi_4_pegs(6, 1, 4, 2, 3)
formatted_moves = format_moves(moves)
for move in formatted_moves:
    print(move)