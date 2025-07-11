def find_filled_groups(q_symbol, m_val):
    """
    Identifies and prints the nonabelian filled groups of order 2*q^m.

    Args:
        q_symbol (str): A string symbol for the odd prime, e.g., 'q'.
        m_val (int): A natural number for the exponent.
    """

    def generate_partitions(n):
        """
        Generates all integer partitions of a number n in non-increasing order.
        """
        def _partitions_helper(target, max_val):
            if target == 0:
                yield []
                return
            if target < 0 or max_val == 0:
                return

            # Partitions that include max_val
            if target >= max_val:
                for p in _partitions_helper(target - max_val, max_val):
                    yield [max_val] + p
            
            # Partitions that do not include max_val
            for p in _partitions_helper(target, max_val - 1):
                yield p
        
        # Using yield from to pass on the generator's values
        yield from _partitions_helper(n, n)

    print(f"The nonabelian filled groups of order 2*({q_symbol})^{m_val} are the generalized dihedral groups Dih(A),")
    print(f"where A is an abelian group of order ({q_symbol})^{m_val}. The possible structures for A are determined")
    print(f"by the partitions of m = {m_val}.\n")

    partitions = list(generate_partitions(m_val))
    num_groups = len(partitions)

    print(f"There are {num_groups} such non-isomorphic groups, corresponding to the {num_groups} partitions of {m_val}:\n")

    for i, p in enumerate(partitions):
        print(f"Group {i+1}: Based on partition {p}")
        
        # Describe the abelian group A
        A_factors = []
        for lambda_i in p:
            A_factors.append(f"C_({q_symbol}^{lambda_i})")
        A_structure = " x ".join(A_factors)
        print(f"  A = {A_structure}")
        
        # Describe the group G = Dih(A)
        print(f"  G = Dih({A_structure})")
        
        # Print the presentation of G
        generators = [f"x_{j+1}" for j in range(len(p))] + ['y']
        print(f"  Presentation: <{', '.join(generators)} | ")
        
        # Relations for the abelian part
        for j, lambda_j in enumerate(p):
            print(f"    x_{j+1}^({q_symbol}^{lambda_j}) = 1,")
        
        if len(p) > 1:
            print(f"    [x_i, x_j] = 1 for 1 <= i < j <= {len(p)},")

        # Relations for the dihedral part
        print(f"    y^2 = 1,")
        for j in range(len(p)):
            terminator = "" if j == len(p) - 1 else ","
            print(f"    y*x_{j+1}*y^(-1) = x_{j+1}^(-1){terminator}")
        print("  >\n")

if __name__ == '__main__':
    # --- Parameters for the problem ---
    # q is an odd prime, represented by a symbol.
    q_symbol = 'q' 
    # m is a natural number.
    # You can change this value to explore other cases.
    m_val = 3
    
    find_filled_groups(q_symbol, m_val)
