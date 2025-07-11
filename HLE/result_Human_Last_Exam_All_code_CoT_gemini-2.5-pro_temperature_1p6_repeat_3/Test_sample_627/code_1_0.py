def calculate_vogel_upper_bound():
    """
    This function calculates an upper bound for the braid index of the
    three-twist knot (6_1) using Vogel's algorithm.
    """
    print("Vogel's algorithm finds an upper bound for the braid index of a knot.")
    print("The formula is: Upper Bound = (Total Number of Critical Points) / 2")
    print("-" * 60)

    # Step 1: Identify critical points on a standard diagram of the three-twist knot.
    # A standard symmetric diagram has 3 outer 'lobes' (radial maxima) and 3
    # inner 'valleys' (radial minima) when viewed from the center.
    num_maxima = 3
    num_minima = 3

    print(f"Step 1: Counting critical points for the three-twist knot diagram.")
    print(f"Number of local distance maxima (lobes): {num_maxima}")
    print(f"Number of local distance minima (valleys): {num_minima}")
    print("-" * 60)

    # Step 2: Calculate the total number of critical points.
    total_critical_points = num_maxima + num_minima
    
    print("Step 2: Calculate the total number of critical points.")
    print(f"Total Critical Points = {num_maxima} + {num_minima}")
    print(f"Total Critical Points = {total_critical_points}")
    print("-" * 60)
    
    # Step 3: Calculate the upper bound for the braid index.
    # The result must be an integer.
    upper_bound = total_critical_points // 2

    print("Step 3: Calculate the upper bound using Vogel's formula.")
    print(f"Upper Bound = {total_critical_points} / 2")
    print(f"Upper Bound = {upper_bound}")
    print("-" * 60)
    
    print(f"The calculated upper bound for the braid index is {upper_bound}.")

calculate_vogel_upper_bound()