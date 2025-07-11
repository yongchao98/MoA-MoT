import math

def get_partitions(n):
    """
    Computes the integer partitions of n.
    For example, for n=4, the partitions are [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]].
    """
    # Using a memoization cache to store results of subproblems
    cache = {}
    def partitions_recursive(target, max_val):
        if (target, max_val) in cache:
            return cache[(target, max_val)]
        
        if target == 0:
            return [[]]
        if target < 0 or max_val == 0:
            return []
        
        result = []
        # We can either include max_val in the partition or not.
        # Case 1: Include max_val (or a smaller number)
        for i in range(min(target, max_val), 0, -1):
            sub_partitions = partitions_recursive(target - i, i)
            for part in sub_partitions:
                result.append([i] + part)
        
        cache[(target, max_val)] = result
        return result

    return partitions_recursive(n, n)

def describe_filled_groups(q, m):
    """
    Describes the nonabelian filled groups of order 2*q**m.
    
    Args:
        q (int): An odd prime number.
        m (int): A natural number.
    """
    if not isinstance(q, int) or not isinstance(m, int) or m < 1:
        print("Error: q must be an integer prime and m must be a natural number (>= 1).")
        return
    if q == 2 or (q > 2 and any(q % i == 0 for i in range(2, int(math.sqrt(q)) + 1))):
        print(f"Error: q={q} is not an odd prime number.")
        return

    print(f"The nonabelian filled groups of order 2*({q}^{m}) = {2 * q**m} are the generalized dihedral groups G = Q \u22CA C\u2082.")
    print("The action of C\u2082 on Q is by inversion (g \u21A6 g\u207B\u00B9).")
    print(f"Q is an abelian group of order {q}^{m}. The possible structures for Q are:")
    
    partitions = get_partitions(m)
    
    for i, part in enumerate(partitions):
        q_factors = []
        for p_val in part:
            if p_val == 1:
                q_factors.append(f"C_{q}")
            else:
                q_factors.append(f"C_({q}^{p_val})")
        
        q_structure = " \u00D7 ".join(q_factors)
        print(f"\n{i+1}. Q \u2245 {q_structure}")
        print(f"   Corresponding filled group G: ({q_structure}) \u22CA C\u2082")

if __name__ == '__main__':
    # Example: Find nonabelian filled groups of order 2 * 3^4 = 162
    q_example = 3
    m_example = 4
    describe_filled_groups(q_example, m_example)
    
    # Another example: order 2 * 5^3 = 250
    # q_example = 5
    # m_example = 3
    # describe_filled_groups(q_example, m_example)
