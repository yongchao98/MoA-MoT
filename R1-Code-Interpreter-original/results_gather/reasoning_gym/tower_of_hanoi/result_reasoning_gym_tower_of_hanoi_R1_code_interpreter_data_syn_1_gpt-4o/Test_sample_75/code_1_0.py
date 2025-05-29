def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    
    moves = []
    k = n - 1  # Optimal number of disks to move to auxiliary peg
    # Move k disks to auxiliary1
    moves += hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
    # Move remaining disk to target using auxiliary2
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    # Move k disks from auxiliary1 to target
    moves += hanoi_4_pegs(k, auxiliary1, target, source, auxiliary2)
    
    return moves

# Solve the problem with 4 disks and 4 pegs
moves = hanoi_4_pegs(4, 4, 2, 1, 3)
print("\n".join(moves))