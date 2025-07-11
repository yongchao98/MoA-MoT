def solve_equilibrium():
    """
    Determines if a body is in equilibrium by summing force vectors.
    """
    try:
        n = int(input())
    except (ValueError, EOFError):
        # Handle cases with no or invalid input for n
        return

    x_components = []
    y_components = []
    z_components = []

    for _ in range(n):
        try:
            x, y, z = map(int, input().split())
            x_components.append(x)
            y_components.append(y)
            z_components.append(z)
        except (ValueError, EOFError):
            # Handle malformed lines or premature end of input
            continue

    sum_x = sum(x_components)
    sum_y = sum(y_components)
    sum_z = sum(z_components)
    
    # Helper to format numbers in the equation string
    def format_component(c):
        return f"({c})" if c < 0 else str(c)

    # Output each number in the final equation as requested
    if x_components:
        print(" + ".join(map(format_component, x_components)), f"= {sum_x}")
        print(" + ".join(map(format_component, y_components)), f"= {sum_y}")
        print(" + ".join(map(format_component, z_components)), f"= {sum_z}")

    # Determine if the body is in equilibrium
    if sum_x == 0 and sum_y == 0 and sum_z == 0:
        print("YES")
    else:
        print("NO")

solve_equilibrium()