def hanoi_4_pegs(n, source, target, auxiliary1, auxiliary2, moves):
    if n == 0:
        return
    if n == 1:
        moves.append(f"Move disk 1 from Peg {source} to Peg {target}")
        return
    hanoi_4_pegs(n - 2, source, auxiliary1, auxiliary2, target, moves)
    moves.append(f"Move disk {n-1} from Peg {source} to Peg {auxiliary2}")
    moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
    moves.append(f"Move disk {n-1} from Peg {auxiliary2} to Peg {target}")
    hanoi_4_pegs(n - 2, auxiliary1, target, source, auxiliary2, moves)

moves = []
hanoi_4_pegs(3, 1, 3, 2, 4, moves)
print("\n".join(moves))