def hanoi(n, source, target, auxiliary):
    if n == 1:
        return [f"Move disk 1 from Peg {source} to Peg {target}"]
    else:
        moves = hanoi(n-1, source, auxiliary, target)
        moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
        moves.extend(hanoi(n-1, auxiliary, target, source))
        return moves

moves = hanoi(6, 3, 2, 1)
print("\n".join(moves))