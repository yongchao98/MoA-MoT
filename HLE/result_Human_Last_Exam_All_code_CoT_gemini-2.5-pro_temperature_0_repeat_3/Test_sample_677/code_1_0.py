def find_max_planar_area():
    """
    Finds the maximum area nm of a chessboard for which the super-knight graph
    is planar, under the assumption that the condition (n-5)(m-5) <= 11 is
    sufficient and n, m >= 6.
    """
    max_area = 0
    best_n, best_m = 0, 0

    # We assume n, m >= 6, so n-5 >= 1 and m-5 >= 1.
    # Let x = n-5 and y = m-5. We need xy <= 11, with x, y >= 1.
    # We can iterate through possible values of n and find the max m for each.
    # The loop for n can stop when (n-5)*(n-5) > 11, as m would have to be < n.
    # (n-5)^2 > 11 => n-5 > sqrt(11) ~= 3.31 => n > 8.31.
    # So we only need to check n = 6, 7, 8.
    
    print("Searching for the maximum planar area nm...")
    print("Condition for planarity: (n-5)(m-5) <= 11, with n, m >= 6\n")

    for n in range(6, 9):
        # For a given n, find the largest m >= 6 satisfying the condition.
        # (m-5) <= 11 / (n-5)
        # m <= 11 / (n-5) + 5
        m = int(11 / (n - 5)) + 5
        
        area = n * m
        print(f"Checking n = {n}:")
        print(f"  Max m such that (n-5)(m-5) <= 11 is m = {m}.")
        print(f"  Area = {n} * {m} = {area}")

        if area > max_area:
            max_area = area
            best_n = n
            best_m = m
    
    # By symmetry, we would get the same results for n > 8.
    # For example, n=16 => 11*(m-5)<=11 => m-5<=1 => m<=6. Area = 16*6=96.
    
    print("\n--- Conclusion ---")
    print(f"The maximum area is found for the rectangle {best_n} x {best_m} (or {best_m} x {best_n}).")
    print("The final equation for the maximum area is:")
    # We want to print the equation with the smaller dimension first.
    if best_n > best_m:
        best_n, best_m = best_m, best_n
    print(f"{best_n} * {best_m} = {max_area}")
    
    return max_area

# Run the function to get the result.
# The function prints the step-by-step reasoning and the final answer.
max_val = find_max_planar_area()

# The final answer is the numerical value.
# print(f"\nFinal Answer: {max_val}")