def hanoi(n, source, auxiliary, destination, moves):
    if n == 1:
        moves.append(f"Move disk 1 from Peg {source} to Peg {destination}")
    else:
        hanoi(n-1, source, destination, auxiliary, moves)
        moves.append(f"Move disk {n} from Peg {source} to Peg {destination}")
        hanoi(n-1, auxiliary, source, destination, moves)

moves = []
hanoi(6, 2, 1, 3, moves)

# Print the moves in chunks to avoid truncation
chunk_size = 20
for i in range(0, len(moves), chunk_size):
    print("\n".join(moves[i:i+chunk_size]))