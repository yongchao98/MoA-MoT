import math

def solve():
    """
    Calculates the number of power subgroups in the generalized quaternion group of size 128 (Q_128).
    
    The group Q_128 has the presentation: <x, y | x^64 = 1, y^2 = x^32, y^-1*x*y = x^-1>.
    A power subgroup G^n is defined as G^n = <g^n | g in G>.
    """

    print("Step 1: Analyze the power subgroups G^n for odd n.")
    print("For any odd integer n, the set of n-th powers {g^n | g in G} generates the entire group G.")
    print("Therefore, G^n = G for all odd n.")
    num_from_odd_n = 1
    print(f"This gives {num_from_odd_n} unique power subgroup (G itself).\n")
    
    print("Step 2: Analyze the power subgroups G^n for even n.")
    print("For any even integer n, the power subgroup G^n is a cyclic subgroup of the form <x^d>.")
    print("Here, d = gcd(n, 64).")
    print("To find the number of distinct subgroups, we need to find all possible values for d.\n")

    order_of_x = 64
    print(f"Since n is even, d = gcd(n, {order_of_x}) must also be an even number that divides {order_of_x}.")
    print(f"So, we need to find all the even divisors of {order_of_x}.")
    
    divisors = [d for d in range(1, order_of_x + 1) if order_of_x % d == 0]
    even_divisors = [d for d in divisors if d % 2 == 0]
    
    print(f"The divisors of {order_of_x} are: {divisors}.")
    print(f"The even divisors of {order_of_x} are: {even_divisors}.")

    num_from_even_n = len(even_divisors)
    print(f"Each of these {num_from_even_n} even divisors corresponds to a unique power subgroup: <x^2>, <x^4>, ..., <x^64>.")
    print(f"This gives {num_from_even_n} unique power subgroups.\n")

    print("Step 3: Calculate the total number of power subgroups.")
    total_subgroups = num_from_odd_n + num_from_even_n
    
    print("The total is the sum of the counts from the odd and even cases.")
    print(f"Total = (subgroups from odd n) + (subgroups from even n)")
    print(f"Total = {num_from_odd_n} + {num_from_even_n} = {total_subgroups}")
    
    return total_subgroups

if __name__ == '__main__':
    solve()
