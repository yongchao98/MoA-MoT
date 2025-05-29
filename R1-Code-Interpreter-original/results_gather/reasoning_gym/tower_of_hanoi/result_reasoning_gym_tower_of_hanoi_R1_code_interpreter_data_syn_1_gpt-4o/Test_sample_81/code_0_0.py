def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    
    moves = []
    # Move n-2 disks from source to auxiliary1 using auxiliary2 and target
    moves += hanoi_4_pegs(n-2, source, auxiliary1, auxiliary2, target)
    # Move the (n-1)th disk to auxiliary2
    moves.append((n-1, source, auxiliary2))
    # Move the nth disk to target
    moves.append((n, source, target))
    # Move the (n-1)th disk from auxiliary2 to target
    moves.append((n-1, auxiliary2, target))
    # Move n-2 disks from auxiliary1 to target using source and auxiliary2
    moves += hanoi_4_pegs(n-2, auxiliary1, target, source, auxiliary2)
    
    return moves

# Generate the sequence of moves for 3 disks and 4 pegs
moves = hanoi_4_pegs(3, 1, 2, 3, 4)
# Format the moves as required
formatted_moves = [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]
print("\n".join(formatted_moves))