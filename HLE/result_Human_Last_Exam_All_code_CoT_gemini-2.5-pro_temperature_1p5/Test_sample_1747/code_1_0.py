def solve_rigid_modules():
    """
    Calculates the number of regular rigid indecomposable modules for a
    complex path algebra of type A_tilde_{2,3}.
    """
    # The type A_tilde_{p,q} corresponds to a canonical algebra with two
    # exceptional tubes of ranks p and q.
    ranks = [2, 3]

    print(f"The algebra of type A_tilde_{ranks[0]},{ranks[1]} has {len(ranks)} exceptional tubes.")
    print("-" * 30)

    # For an exceptional tube of rank r, there are r-1 regular rigid
    # indecomposable modules.
    num_rigid_modules_per_tube = []

    # Calculate for the first tube
    rank1 = ranks[0]
    rigid1 = rank1 - 1
    num_rigid_modules_per_tube.append(rigid1)
    print(f"The rank of the first tube is {rank1}.")
    print(f"Number of regular rigid indecomposable modules in this tube: {rank1} - 1 = {rigid1}")
    print("-" * 30)

    # Calculate for the second tube
    rank2 = ranks[1]
    rigid2 = rank2 - 1
    num_rigid_modules_per_tube.append(rigid2)
    print(f"The rank of the second tube is {rank2}.")
    print(f"Number of regular rigid indecomposable modules in this tube: {rank2} - 1 = {rigid2}")
    print("-" * 30)

    # Calculate the total
    total_rigid_modules = sum(num_rigid_modules_per_tube)
    
    # Print the final equation with all numbers
    num1 = num_rigid_modules_per_tube[0]
    num2 = num_rigid_modules_per_tube[1]
    total = total_rigid_modules
    print(f"Total number of regular rigid indecomposable modules: {num1} + {num2} = {total}")

solve_rigid_modules()