import math

def get_divisors(n):
    """Returns a list of all divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return list(divs)

def phi(n):
    """Calculates Euler's totient function phi(n)."""
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

def solve_graph_counting():
    """
    Calculates the number of non-isomorphic, connected, 3-regular, adjustable graphs
    with 2000 vertices that have a perfect matching.
    This is equivalent to counting binary necklaces of length 1000.
    """
    n = 1000

    # The number of non-isomorphic graphs is calculated using Burnside's Lemma for the
    # action of the dihedral group D_n on binary strings of length n.
    # The number of orbits is (1/|D_n|) * sum_{g in D_n} |X^g|,
    # where X is the set of 2^n strings and |X^g| is the number of strings fixed by g.
    
    # 1. Contribution from rotations
    # The sum of fixed points for all rotations is given by: sum_{d|n} phi(n/d) * 2^d
    divisors = get_divisors(n)
    rotation_sum = 0
    for d in divisors:
        rotation_sum += phi(n // d) * (2**d)

    # 2. Contribution from reflections
    # For n=1000 (even), let m = n/2 = 500.
    m = n // 2
    # There are m reflections through opposite vertices. Each fixes 2^(m+1) strings.
    reflection_sum_1 = m * (2**(m + 1))
    # There are m reflections through midpoints of opposite edges. Each fixes 2^m strings.
    reflection_sum_2 = m * (2**m)
    reflection_sum = reflection_sum_1 + reflection_sum_2

    # 3. Total calculation
    # |D_n| = 2*n
    group_size = 2 * n
    total_fixed_configs = rotation_sum + reflection_sum
    num_graphs = total_fixed_configs // group_size

    # The problem asks to output the numbers in the final equation.
    print("The number of such graphs is found by counting binary necklaces of length 1000.")
    print("Using Burnside's Lemma, the formula is (sum of fixed points) / (size of group).")
    print("\nCalculation steps:")
    print(f"1. Sum of configurations fixed by rotations: {rotation_sum}")
    print(f"2. Sum of configurations fixed by reflections: {reflection_sum}")
    print(f"3. Total fixed configurations: {rotation_sum} + {reflection_sum} = {total_fixed_configs}")
    print(f"4. Size of the dihedral group D_{n}: {group_size}")

    print("\nThe final equation is:")
    print(f"({rotation_sum} + {reflection_sum}) / {group_size} = {num_graphs}")
    
    return num_graphs

if __name__ == '__main__':
    final_answer = solve_graph_counting()
    # The final answer is an extremely large number.
    # We will wrap it in the required format.
    # print(f"\n<<< {final_answer} >>>")

solve_graph_counting()