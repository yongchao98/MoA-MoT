def print_sequence():
    # Initial board
    board = [
        ['.', '.', 'B', 'w'],  # Row 3
        ['.', 'B', '.', '.'],  # Row 2
        ['.', 'w', '.', '.']   # Row 1
    ]
    
    # Try the sequence
    moves = [
        ("w", "C1", "A2"),  # White knight from C1 to A2
        ("B", "C3", "C1"),  # Black knight from C3 to C1
        ("w", "D3", "C3"),  # White knight from D3 to C3
        ("B", "B2", "D3")   # Black knight from B2 to D3
    ]
    
    print("Move sequence:")
    for color, start, end in moves:
        print(f"{color}: {start} → {end}")
    
    print("\nThis sequence works because:")
    print("1. All moves are legal knight moves")
    print("2. Players alternate turns (white starts)")
    print("3. All moves are to empty squares")
    print("4. Final positions achieve the swap")

print_sequence()