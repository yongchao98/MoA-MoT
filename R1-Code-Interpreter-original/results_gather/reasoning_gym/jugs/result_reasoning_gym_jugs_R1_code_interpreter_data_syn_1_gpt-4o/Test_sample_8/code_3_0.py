moves = [
    "fill A",       # Step 1: Fill Jug A (11 liters)
    "pour A->B",    # Step 2: Pour from A to B (Jug A: 5 liters, Jug B: 6 liters)
    "empty B",      # Step 3: Empty Jug B
    "pour A->B",    # Step 4: Pour from A to B (Jug A: 0 liters, Jug B: 5 liters)
    "fill A",       # Step 5: Fill Jug A (11 liters)
    "pour A->B"     # Step 6: Pour from A to B (Jug A: 8 liters, Jug B: 6 liters)
]

print(moves)