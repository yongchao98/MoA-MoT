def describe_sets_for_n_components(n):
    """
    This function describes two closed, connected sets A and B whose union is the
    unit square [0,1]x[0,1] and whose intersection has n connected components.

    Args:
        n (int): The desired number of components in the intersection.
                 Must be a positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        print("Please provide a positive integer for the number of components.")
        return

    print(f"Constructing sets for n = {n} intersection components.\n")

    # The height of each horizontal stripe
    h = 1.0 / (n + 1)

    # Describe Set A
    print("Set A is the union of:")
    print("1. The vertical line segment {(0, y) | 0 <= y <= 1}.")
    a_stripes = []
    # A takes the stripes starting from k=1 (the second stripe from the bottom)
    for k in range(n + 1):
        if k % 2 == 1:
            a_stripes.append(f"{(k)*h:.4f} <= y <= {(k+1)*h:.4f}")
    print(f"2. The horizontal stripes {{ (x, y) | 0 <= x <= 1 and y is in one of the intervals: }}")
    for s in a_stripes:
        print(f"   - [y: {s}]")
    print("\nSet A is connected and closed.\n")


    # Describe Set B
    print("Set B is the union of:")
    print("1. The vertical line segment {(1, y) | 0 <= y <= 1}.")
    b_stripes = []
    # B takes the stripes starting from k=0 (the bottom stripe)
    for k in range(n + 1):
        if k % 2 == 0:
            b_stripes.append(f"{(k)*h:.4f} <= y <= {(k+1)*h:.4f}")
    print(f"2. The horizontal stripes {{ (x, y) | 0 <= x <= 1 and y is in one of the intervals: }}")
    for s in b_stripes:
        print(f"   - [y: {s}]")
    print("\nSet B is connected and closed.\n")

    # Describe the intersection A ∩ B
    print("-" * 20)
    print("The union A U B is the unit square.")
    print(f"The intersection A ∩ B consists of {n} connected components.")
    print("These components are the horizontal lines separating the stripes:\n")
    
    # The "final equation" is the list of components
    for i in range(1, n + 1):
        y_val = i * h
        print(f"Component {i}: The line segment {{ (x, {y_val:.4f}) | 0 <= x <= 1 }}")

# You can change this value to see the construction for any number of components.
number_of_components = 5
describe_sets_for_n_components(number_of_components)