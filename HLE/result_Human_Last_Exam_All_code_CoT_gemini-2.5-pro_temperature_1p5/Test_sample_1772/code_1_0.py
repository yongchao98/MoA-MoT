def generate_sets(k=10):
    """
    Generates and prints the first k elements of the example sets A and B.
    """
    # Set A = {0} U {1/n | n >= 1}
    set_A = {0.0}
    for n in range(1, k + 1):
        set_A.add(1/n)

    # Set B = A U {m | m >= 2}
    set_B = set_A.copy()
    for m in range(2, k + 2):
        set_B.add(float(m))

    print("Sample of elements from Set A:")
    # Sort for clear presentation
    print(sorted(list(set_A)))

    print("\nSample of elements from Set B:")
    print(sorted(list(set_B)))

generate_sets()