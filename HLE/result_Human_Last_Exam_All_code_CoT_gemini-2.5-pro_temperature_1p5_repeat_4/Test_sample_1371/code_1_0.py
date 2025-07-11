import math

def solve_seating_arrangement():
    """
    Calculates the number of ways to seat the attendees at a circular table
    based on the given constraints.
    """

    # Step 1: Calculate the number of external arrangements.
    # We are arranging 5 items in a circle: Block A (Academics), Block C' (4 Classicists),
    # and 3 'buffer' individuals (2 Ethicists, 1 Cassie).
    # The constraint is that Block A cannot be adjacent to Block C'.
    # Total arrangements of 5 items in a circle = (5-1)! = 4!
    # Forbidden arrangements (A next to C'): Treat (AC') as one item.
    # Arrangements of 4 items in a circle = (4-1)! = 3!. A and C' can swap, so *2.
    total_items = 5
    num_external = math.factorial(total_items - 1) - (math.factorial(total_items - 2) * 2)

    # Step 2: Calculate the internal permutations for the Academic Block (A).
    # Block A consists of 12 scientists and 4 mathematicians.
    # It can be ordered as (Scientists - Mathematicians) or (Mathematicians - Scientists), hence factor of 2.
    # Scientists: 10 non-rowers (10!) and 2 rowers (2!).
    # Mathematicians: 3 non-rowers (3!).
    order_factor = 2
    perms_non_rower_scientists = math.factorial(10)
    perms_rower_scientists = math.factorial(2)
    perms_non_rower_mathematicians = math.factorial(3)
    num_internal_A = order_factor * perms_non_rower_scientists * perms_rower_scientists * perms_non_rower_mathematicians

    # Step 3: Calculate the internal permutations for the Classicist Block (C').
    # This block has 4 non-Cassie classicists.
    num_internal_C = math.factorial(4)

    # Step 4: Calculate the total number of arrangements.
    total_arrangements = num_external * num_internal_A * num_internal_C
    
    # Print the breakdown of the calculation
    print("The final calculation is based on multiplying the number of ways to arrange the groups by the number of ways to arrange the people within each group.")
    print(f"\n1. Number of ways to arrange the blocks and buffer individuals externally: {num_external}")
    
    print(f"\n2. Number of ways to internally arrange the Academic Block (A):")
    print(f"   ({order_factor} orders) * ({perms_non_rower_scientists} for non-rower scientists) * ({perms_rower_scientists} for rower scientists) * ({perms_non_rower_mathematicians} for non-rower mathematicians) = {num_internal_A}")

    print(f"\n3. Number of ways to internally arrange the Classicist Block (C'): {num_internal_C}")

    print("\nTotal number of ways to arrange the table is:")
    print(f"{num_external} (External) * {num_internal_A} (Academic Block) * {num_internal_C} (Classicist Block) = {total_arrangements:,}")


solve_seating_arrangement()
<<<25082265600>>>