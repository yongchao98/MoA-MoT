def hanoi(n, source, target, auxiliary, moves):
    if n == 1:
        moves.append(f"Move disk 1 from Peg {source} to Peg {target}")
    else:
        hanoi(n-1, source, auxiliary, target, moves)
        moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
        hanoi(n-1, auxiliary, target, source, moves)

moves = []
hanoi(6, 3, 1, 2, moves)

# Print the moves in chunks to avoid truncation
for i in range(0, len(moves), 20):
    print("\n".join(moves[i:i+20]))
    print("---")