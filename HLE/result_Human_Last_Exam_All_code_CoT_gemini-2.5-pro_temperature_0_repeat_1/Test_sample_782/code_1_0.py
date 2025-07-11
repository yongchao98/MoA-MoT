def calculate_cost():
    """
    Calculates the minimum cost for computing 2A - 3B on a twisted Edwards curve.

    The calculation is based on the identity 2A - 3B = 2(A - B) - B.
    - A, B are in affine coordinates.
    - The result is in extended coordinates.
    - Cost is in field multiplications (M), with S=M.
    """

    # Step 1: Cost of C = A - B
    # This involves adding two affine points (A and -B) to get an extended point.
    # Cost = 2M (for T coordinates) + 8M (for addition)
    cost_step1 = 10

    # Step 2: Cost of D = 2C
    # This is a doubling of a general extended point.
    # Cost = 4S + 5M = 9M
    cost_step2 = 9

    # Step 3: Cost of R = D - B
    # This is a mixed addition of an extended point and an affine point.
    # Cost = 10M
    cost_step3 = 10

    # Total cost is the sum of the costs of each step.
    total_cost = cost_step1 + cost_step2 + cost_step3

    print("The calculation is performed using the identity: 2A - 3B = 2(A - B) - B")
    print("-" * 60)
    print(f"Cost of Step 1 (computing A - B): {cost_step1}M")
    print(f"Cost of Step 2 (doubling the result): {cost_step2}M")
    print(f"Cost of Step 3 (subtracting B): {cost_step3}M")
    print("-" * 60)
    print("The final cost equation is:")
    print(f"Total Cost = {cost_step1} + {cost_step2} + {cost_step3} = {total_cost} multiplications")

calculate_cost()