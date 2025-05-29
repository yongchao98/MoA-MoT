def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    moves = []
    if n == 2:
        # Directly solve for 2 disks
        moves.append((1, source, auxiliary1))
        moves.append((2, source, target))
        moves.append((1, auxiliary1, target))
    else:
        k = 2  # A choice for k, can be optimized further
        # Move k disks to auxiliary1
        moves += hanoi_4_pegs(k, source, auxiliary1, auxiliary2, target)
        # Move remaining n-k disks to target using auxiliary2
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
    # Move n-1 disks to auxiliary
    moves += hanoi_3_pegs(n-1, source, auxiliary, target)
    # Move nth disk to target
    moves.append((n, source, target))
    # Move n-1 disks from auxiliary to target
    moves += hanoi_3_pegs(n-1, auxiliary, target, source)
    return moves

# Generate the sequence of moves for 4 disks and 4 pegs
moves = hanoi_4_pegs(4, 3, 1, 2, 4)

# Format the output
formatted_moves = [f"Move disk {disk} from Peg {src} to Peg {tgt}" for disk, src, tgt in moves]
print("\n".join(formatted_moves))