import math

def combinations(n, k):
    """
    Calculates the number of combinations of choosing k items from a set of n.
    This is also known as "n choose k".
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    # Optimization: C(n, k) == C(n, n-k)
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        # Efficiently calculate n! / (k! * (n-k)!)
        res = res * (n - i) // (i + 1)
    return res

def solve_pizza_combinations():
    """
    Calculates the number of possible pizza size combinations based on the user's constraints.
    """
    pizzas_per_area = {}

    # 1. Iterate through all possible rounded slice areas (from 5.0 to 20.0)
    # Using integers for the loop (50 to 200) avoids floating-point inaccuracies.
    for target_a_int in range(50, 201):
        target_a = target_a_int / 10.0
        
        valid_diameters_for_area = []
        
        # 2. Iterate through all possible integer diameters (1 to 20 inches)
        for d in range(1, 21):
            pizza_area = math.pi * (d / 2)**2
            
            # Find the range for the number of slices 'N' so that round(pizza_area / N, 1) == target_a.
            # This is equivalent to: target_a - 0.05 <= pizza_area / N < target_a + 0.05
            # We rearrange to solve for N:
            # pizza_area / (target_a + 0.05) < N <= pizza_area / (target_a - 0.05)
            
            # Calculate lower and upper bounds for N
            lower_bound_n = pizza_area / (target_a + 0.05)
            upper_bound_n = pizza_area / (target_a - 0.05) if target_a > 5.0 else float('inf')

            # Find the smallest integer N strictly greater than the lower bound.
            min_n = math.floor(lower_bound_n) + 1
            
            # The number of slices N must be a multiple of 4. Find the first valid candidate for N.
            if min_n < 4:
                first_valid_n = 4
            else:
                rem = min_n % 4
                if rem == 0:
                    first_valid_n = min_n
                else:
                    first_valid_n = min_n + (4 - rem)
            
            # If a valid N exists for this diameter 'd', add 'd' to our list.
            if first_valid_n <= upper_bound_n:
                valid_diameters_for_area.append(d)

        # 3. Store the list of valid diameters for the current target area.
        if valid_diameters_for_area:
            pizzas_per_area[target_a] = valid_diameters_for_area

    # 4. Calculate the number of combinations for sets of 3, 4, and 5 pizzas.
    count3 = 0
    count4 = 0
    count5 = 0

    for area in pizzas_per_area:
        k = len(pizzas_per_area[area])
        if k >= 3:
            count3 += combinations(k, 3)
        if k >= 4:
            count4 += combinations(k, 4)
        if k >= 5:
            count5 += combinations(k, 5)
            
    total_count = count3 + count4 + count5
    
    # 5. Print the final results.
    print(f"Number of possible sets with 3 pizza sizes: {count3}")
    print(f"Number of possible sets with 4 pizza sizes: {count4}")
    print(f"Number of possible sets with 5 pizza sizes: {count5}")
    print(f"Total number of possible sets: {count3} + {count4} + {count5} = {total_count}")

# Execute the main function to find the answer.
solve_pizza_combinations()
<<<2740>>>