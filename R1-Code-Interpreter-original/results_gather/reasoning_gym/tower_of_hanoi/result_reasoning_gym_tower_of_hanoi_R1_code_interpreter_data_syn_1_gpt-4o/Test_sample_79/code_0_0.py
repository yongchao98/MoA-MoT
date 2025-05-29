def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = []
    k = n - 2  # Optimal number of disks to move to auxiliary peg
    # Move top k disks to auxiliary1
    moves += hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
    # Move remaining n-k disks to target
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
    # Move top n-1 disks to auxiliary
    moves += hanoi_3_pegs(n-1, source, auxiliary, target)
    # Move nth disk to target
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    # Move n-1 disks from auxiliary to target
    moves += hanoi_3_pegs(n-1, auxiliary, target, source)
    
    return moves

# Generate the sequence of moves for 4 disks and 4 pegs
moves = hanoi_4_pegs(4, 4, 1, 2, 3)
print("\n".join(moves))