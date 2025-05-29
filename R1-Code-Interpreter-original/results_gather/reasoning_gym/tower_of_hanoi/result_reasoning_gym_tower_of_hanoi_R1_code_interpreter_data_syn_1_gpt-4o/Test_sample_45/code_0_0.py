def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    # Find the optimal number of disks to move to the auxiliary peg
    k = n - 2
    moves = []
    
    # Move k disks to auxiliary1
    moves += hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
    
    # Move the remaining n-k disks to target using 3 pegs
    moves += hanoi_3_pegs(n-k, source, target, auxiliary2)
    
    # Move k disks from auxiliary1 to target
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    
    return moves

def hanoi_3_pegs(n, source, target, auxiliary):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = []
    moves += hanoi_3_pegs(n-1, source, auxiliary, target)
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    moves += hanoi_3_pegs(n-1, auxiliary, target, source)
    
    return moves

# Solve the problem for 5 disks, moving from Peg 2 to Peg 1
moves = hanoi_4_pegs(5, 2, 1, 3, 4)
for move in moves:
    print(move)