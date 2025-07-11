def solve_equilibrium():
    """
    Parses the vector input and prints the summation equations
    to show why the body is in equilibrium.
    """
    input_data = """3
3 -1 7
-5 2 -4
2 -1 -3"""

    lines = input_data.strip().split('\n')
    num_vectors = int(lines[0])
    vectors = []
    for i in range(1, num_vectors + 1):
        try:
            parts = list(map(int, lines[i].split()))
            if len(parts) == 3:
                vectors.append(parts)
        except (ValueError, IndexError):
            # Skip malformed lines if any
            continue

    if not vectors:
        print("No valid vectors found in input.")
        return

    x_coords = [v[0] for v in vectors]
    y_coords = [v[1] for v in vectors]
    z_coords = [v[2] for v in vectors]

    def format_and_print_equation(name, coords):
        """Formats and prints a single summation equation."""
        total = sum(coords)
        # Formats numbers like (3) + (-5) + (2) for clarity
        parts = [f"({num})" for num in coords]
        equation_str = " + ".join(parts)
        print(f"Sum of {name}-components: {equation_str} = {total}")

    format_and_print_equation("x", x_coords)
    format_and_print_equation("y", y_coords)
    format_and_print_equation("z", z_coords)

    # Final result check
    if sum(x_coords) == 0 and sum(y_coords) == 0 and sum(z_coords) == 0:
        print("\nResult: The body is in equilibrium (YES).")
    else:
        print("\nResult: The body is not in equilibrium (NO).")

solve_equilibrium()