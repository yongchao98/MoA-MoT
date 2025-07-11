import math

def solve_power_subgroups():
    """
    Calculates the number of power subgroups in the semidihedral group of a given size.
    """
    group_size = 512

    # A semidihedral group SD_{2^n} has order 2^n.
    # We find n for the given group size.
    if group_size <= 0 or (group_size & (group_size - 1)) != 0:
        print("Error: Group size must be a power of 2.")
        return

    n = int(math.log2(group_size))

    if n < 3:
        print(f"Error: The semidihedral group SD_{group_size} is not defined for n < 3.")
        return

    # The exponent of the semidihedral group SD_{2^n} is 2^(n-1).
    exponent_power = n - 1
    exponent = 2**exponent_power

    # The distinct power subgroups are G^d, where d is a divisor of the exponent.
    # The divisors of the exponent 2^(n-1) are 2^0, 2^1, ..., 2^(n-1).
    # The number of these divisors is (n-1) - 0 + 1 = n.
    num_power_subgroups = n

    print(f"For the semidihedral group of size {group_size}, we have n = {n}.")
    print(f"The exponent of this group is 2^({n} - 1) = {exponent}.")
    print(f"The number of power subgroups is the number of divisors of the exponent {exponent}.")
    
    divisors_powers = list(range(n))
    divisors = [2**i for i in divisors_powers]
    
    # We are asked to output each number in the final equation.
    # The number of divisors of 2^(n-1) is (n-1) + 1 = n.
    print(f"Number of divisors = ({n} - 1) + 1 = {n}")
    print(f"The divisors are: {', '.join(map(str, divisors))}.")
    print(f"This leads to {len(divisors)} distinct power subgroups.")
    print(f"Thus, the number of power subgroups is {num_power_subgroups}.")

solve_power_subgroups()