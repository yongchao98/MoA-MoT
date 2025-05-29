def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2):
    if n == 0:
        return []
    if n == 1:
        return [(1, source, target)]
    if n == 2:
        return [(1, source, auxiliary1), (2, source, target), (1, auxiliary1, target)]
    
    moves = []
    moves += hanoi_4_pegs(n-2, source, auxiliary1, auxiliary2, target)
    moves.append((n-1, source, auxiliary2))
    moves.append((n, source, target))
    moves.append((n-1, auxiliary2, target))
    moves += hanoi_4_pegs(n-2, auxiliary1, target, source, auxiliary2)
    return moves

# Generate the sequence of moves for 3 disks and 4 pegs
moves = hanoi_4_pegs(3, 1, 2, 3, 4)

# Format the output
formatted_moves = [f"Move disk {disk} from Peg {src} to Peg {dst}" for disk, src, dst in moves]
print("\n".join(formatted_moves))