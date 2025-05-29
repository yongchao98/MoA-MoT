def hanoi(n, source, target, auxiliary, moves):
    if n == 1:
        moves.append(f"Move disk 1 from Peg {source} to Peg {target}")
    else:
        hanoi(n-1, source, auxiliary, target, moves)
        moves.append(f"Move disk {n} from Peg {source} to Peg {target}")
        hanoi(n-1, auxiliary, target, source, moves)

moves = []
hanoi(6, 2, 1, 3, moves)

# Store the output in a single string
output = "\n".join(moves)

# Print the output in smaller chunks
chunk_size = 1000
for i in range(0, len(output), chunk_size):
    print(output[i:i+chunk_size])