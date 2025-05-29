def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    # Move the top 2 disks to auxiliary1
    moves = hanoi_3_pegs(2, source, auxiliary1, auxiliary2)
    
    # Move the remaining 3 disks to target using 3 pegs
    moves += hanoi_3_pegs(3, source, target, auxiliary2)
    
    # Move the 2 disks from auxiliary1 to target
    moves += hanoi_3_pegs(2, auxiliary1, target, source)
    
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