def hanoi(n, source, target, auxiliary, moves):
    if n == 1:
        moves.append(f"Move disk 1 from Peg {source} to Peg {target}")
    else:
        hanoi(n-1, source, auxiliary, target, moves)
        moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
        hanoi(n-1, auxiliary, target, source, moves)

moves = []
hanoi(4, 2, 1, 3, moves)
print("\n".join(moves))