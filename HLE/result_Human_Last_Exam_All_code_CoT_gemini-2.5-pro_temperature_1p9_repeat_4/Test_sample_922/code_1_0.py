def find_perfect_nth_root(n, p):
    """
    Finds the integer p-th root of n using binary search.
    Returns the integer root if it's a perfect p-th power, otherwise None.
    """
    # Establish a search range. An upper bound of n is safe for n > 1.
    low = 1
    high = n

    while low <= high:
        mid = (low + high) // 2
        if mid == 0:
            low = 1
            continue

        try:
            val = mid ** p
        except OverflowError:
            # If mid**p is too large to compute, it's definitely greater than n.
            val = float('inf')

        if val == n:
            return mid  # Found the perfect integer root
        elif val < n:
            low = mid + 1
        else:  # val > n
            high = mid - 1
            
    return None # No perfect integer root was found

def solve_sequence():
    """
    Solves the sixth powers problem presented in the sequence.
    """
    # The sequence of numbers are the bases for the sixth powers.
    bases = [24663, 35005, 119261, 196219, 211770, 227296]
    power = 6

    # 1. Calculate the sum of the 6th powers.
    sum_of_powers = sum(n ** power for n in bases)

    # 2. Find the integer 6th root of the sum.
    result = find_perfect_nth_root(sum_of_powers, power)

    # 3. If a solution is found, print the full equation.
    if result is not None:
        print("The given sequence consists of the bases in a sum of six sixth powers that equals another sixth power.")
        print("The full equation is:")
        
        # Build and print the equation string, including the final result.
        equation_parts = [f"{n}^{power}" for n in bases]
        equation = " + ".join(equation_parts) + f" = {result}^{power}"
        print(equation)
        
        print(f"\nThe integer that completes the sequence is {result}.")
        # This will be used for the final answer format.
        print(f"\n<<<{result}>>>")
    else:
        print("No integer solution was found.")

solve_sequence()