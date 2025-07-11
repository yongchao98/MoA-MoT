import math

def calculate_and_print_orders():
    """
    This function calculates the minimal order u_r(n) for n from 3 to 12
    and prints each calculation, fulfilling the prompt's request to show
    the numbers in the "final equation" used for the computation.
    """
    print("Calculating the values for u_r(n) from n=3 to 12:")
    
    results_list = []
    for n in range(3, 13):
        if n % 2 != 0:
            # For odd n, the formula is n-1
            result = n - 1
            print(f"For n={n} (odd): u_r({n}) = {n} - 1 = {result}")
            results_list.append(result)
        else:
            # For even n, the formula is 2 * floor(n / 4)
            floor_val = n // 4
            result = 2 * floor_val
            print(f"For n={n} (even): u_r({n}) = 2 * floor({n}/4) = 2 * {floor_val} = {result}")
            results_list.append(result)

    print("\nThe final set of values {u_r(3), u_r(4), ..., u_r(12)} is:")
    print(results_list)

calculate_and_print_orders()