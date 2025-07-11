import math

def solve_power_subgroups():
    """
    Calculates the number of power subgroups in the generalized quaternion group Q_128.
    """
    group_size = 128

    # Step 1: Identify the parameter n for the group Q_{2^n}
    # group_size = 2^n, so n = log2(group_size)
    if group_size <= 0 or (group_size & (group_size - 1)) != 0:
        print("Error: Group size must be a power of 2.")
        return
        
    n = int(math.log2(group_size))
    print(f"The generalized quaternion group of size {group_size} is denoted as Q_{{{group_size}}}, which is Q_{{2^{n}}}.")

    # Step 2: Determine the exponent of the group Q_{2^n}
    # The exponent for Q_{2^n} (n>=3) is 2^(n-1).
    exponent = 2**(n - 1)
    print(f"The exponent of this group is 2^({n}-1) = {exponent}.")

    # Step 3: Find all divisors of the exponent.
    # The number of distinct power subgroups is equal to the number of divisors of the exponent.
    print(f"The number of power subgroups is equal to the number of divisors of the exponent, {exponent}.")
    
    divisors = []
    for i in range(1, exponent + 1):
        if exponent % i == 0:
            divisors.append(i)
    
    num_divisors = len(divisors)
    
    print(f"The divisors of {exponent} are: {', '.join(map(str, divisors))}.")

    # Step 4: Display the final count as an equation.
    # The final equation shows a '1' for each divisor found, summing to the total count.
    equation_str = " + ".join(['1'] * num_divisors)
    print("Each divisor corresponds to a unique power subgroup. Counting them gives the equation:")
    print(f"{equation_str} = {num_divisors}")
    
    print(f"\nThus, there are {num_divisors} power subgroups in Q_{group_size}.")

solve_power_subgroups()