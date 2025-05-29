def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = []
    k = 2  # Optimal number of disks to move to an intermediate peg for n=5
    # Move k disks from source to auxiliary1 using all pegs
    moves += hanoi_4_pegs(k, source, auxiliary1, target, auxiliary2)
    # Move remaining disks from source to target using auxiliary2
    moves += hanoi_3_pegs(n - k, source, target, auxiliary2)
    # Move k disks from auxiliary1 to target using all pegs
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    return moves

def hanoi_3_pegs(n, source, target, auxiliary):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = []
    moves += hanoi_3_pegs(n - 1, source, auxiliary, target)
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    moves += hanoi_3_pegs(n - 1, auxiliary, target, source)
    return moves

# Solve the problem with 5 disks, moving from Peg 2 to Peg 4
moves = hanoi_4_pegs(5, 2, 4, 1, 3)
for move in moves:
    print(move)