def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    moves = []
    # Optimal choice of m for 6 disks and 4 pegs is known to be 3
    m = 3
    # Move m disks to auxiliary1 using all pegs
    moves += hanoi_4_pegs(m, source, auxiliary1, auxiliary2, target)
    # Move the remaining n-m disks directly to target using auxiliary2
    moves += hanoi_3_pegs(n - m, source, target, auxiliary2)
    # Move m disks from auxiliary1 to target using all pegs
    moves += hanoi_4_pegs(m, auxiliary1, target, source, auxiliary2)
    
    return moves

def hanoi_3_pegs(n, source, target, auxiliary):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    moves = []
    # Move n-1 disks to auxiliary
    moves += hanoi_3_pegs(n - 1, source, auxiliary, target)
    # Move the nth disk to target
    moves.append((n, source, target))
    # Move n-1 disks from auxiliary to target
    moves += hanoi_3_pegs(n - 1, auxiliary, target, source)
    
    return moves

# Generate the sequence of moves for 6 disks from Peg 2 to Peg 1
moves = hanoi_4_pegs(6, 2, 1, 3, 4)

# Format the output
formatted_moves = [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]
print("\n".join(formatted_moves))