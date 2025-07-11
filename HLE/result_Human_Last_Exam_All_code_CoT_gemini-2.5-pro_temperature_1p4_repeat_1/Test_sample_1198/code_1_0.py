import math

def check_construction_mod3():
    """
    This function attempts to construct an interval for a number 'a' such that
    floor(a^n) == n (mod 3) for all n > 0.
    It starts with the simplest case: floor(a) = 1.
    """
    # Initialize the interval for 'a' based on floor(a) = 1.
    # We use a list [low, high] to represent the interval [low, high).
    # k_n represents floor(a^n)
    k1 = 1
    # Check if k1 satisfies the condition for n=1
    if k1 % 3 != 1:
        print(f"Initial choice k1={k1} is invalid.")
        return

    # a is in the interval [k1, k1 + 1)
    a_interval = [float(k1), float(k1 + 1)]
    print(f"n=1: Let floor(a) = {k1}. This implies a is in [{a_interval[0]}, {a_interval[1]}).")

    for n in range(2, 10):
        # The condition for the n-th step
        target_rem = n % 3
        
        # Calculate the range for a^n based on the current interval for 'a'
        # low_bound <= a^n < high_bound
        low_bound = a_interval[0] ** n
        high_bound = a_interval[1] ** n
        
        # Determine the possible integer values for floor(a^n)
        # Small tolerance for floating point inaccuracies
        epsilon = 1e-9
        min_k_n = math.ceil(low_bound - epsilon)
        max_k_n = math.floor(high_bound - epsilon)
        
        possible_k_n = []
        for k in range(min_k_n, max_k_n + 1):
            if k % 3 == target_rem:
                possible_k_n.append(k)
        
        print(f"\nn={n}: Need floor(a^n) === {target_rem} (mod 3).")
        print(f"Current interval for 'a' is [{a_interval[0]:.4f}, {a_interval[1]:.4f}).")
        print(f"This implies a^n is in [{low_bound:.4f}, {high_bound:.4f}).")
        print(f"Possible integer values for floor(a^n) are from {min_k_n} to {max_k_n}.")

        if not possible_k_n:
            print(f"\nFAILURE at n={n}:")
            print(f"None of the possible integers for floor(a^n) (from {min_k_n} to {max_k_n}) satisfy the condition of being congruent to {target_rem} modulo 3.")
            print("This construction path has failed.")
            return
            
        print(f"Integers in this range satisfying the modulo 3 condition: {possible_k_n}")
        
        # For this demonstration, we pick the first possible k_n.
        # A full proof would need to show that all branches fail.
        # However, in the first few steps, the choice is forced.
        if len(possible_k_n) > 1:
            print(f"Multiple choices available for floor(a^n). We will proceed with the smallest one: {possible_k_n[0]}.")

        k_n = possible_k_n[0]

        # Update the interval for 'a'
        # a is in [k_n^(1/n), (k_n+1)^(1/n))
        new_low = k_n**(1/n)
        new_high = (k_n+1)**(1/n)
        
        a_interval[0] = max(a_interval[0], new_low)
        a_interval[1] = min(a_interval[1], new_high)

        if a_interval[0] >= a_interval[1]:
            print(f"FAILURE at n={n}: The interval for 'a' became empty.")
            return

    print("\nConstruction succeeded for the tested range.")

check_construction_mod3()