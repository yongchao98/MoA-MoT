import math

def get_integer_partitions(n):
    """
    Generates all integer partitions of n.
    Credit: Adapted from a recursive algorithm by Jerome Kelleher.
    """
    partitions = []
    def partition_recursive(target, max_val, current_partition):
        if target == 0:
            partitions.append(current_partition)
            return
        if max_val == 0:
            return

        for i in range(min(max_val, target), 0, -1):
            partition_recursive(target - i, i, current_partition + [i])

    partition_recursive(n, n, [])
    return partitions

def describe_group(q, partition):
    """
    Describes the abelian group Q based on a partition.
    Example: q=3, partition=[2, 1] -> C_9 x C_3
    """
    if not partition:
        return "C_1" # Trivial group
    
    parts = [f"C_{q**p}" for p in partition]
    return " x ".join(parts)

def find_filled_groups(q, m):
    """
    Finds and describes the nonabelian filled groups of order 2*q**m.
    
    Args:
        q (int): An odd prime number.
        m (int): A natural number.
    """
    
    if not (isinstance(q, int) and q > 2 and all(q % i != 0 for i in range(2, int(math.sqrt(q)) + 1))):
        print(f"Error: q must be an odd prime. Received q={q}.")
        return

    if not (isinstance(m, int) and m >= 1):
        print(f"Error: m must be a natural number (>= 1). Received m={m}.")
        return

    print(f"The nonabelian filled groups of order 2*{q}^{m} for q={q}, m={m} are the generalized dihedral groups Dih(Q), where Q is an abelian group of order {q**m}.")
    print("The possible structures for Q are determined by the integer partitions of m:")
    
    partitions = get_integer_partitions(m)
    
    for i, p in enumerate(partitions):
        q_description = describe_group(q, p)
        # The final output format requires printing numbers in the final equation.
        # While the groups are symbolic, we print the order calculation.
        group_order_str = f"2 * {q}^{m} = {2 * q**m}"
        print(f"  {i+1}. Q = {q_description}. Group: Dih({q_description}). Order: {group_order_str}")


# --- Example Usage ---
# You can change these values to explore other cases.
q_example = 3
m_example = 4

find_filled_groups(q_example, m_example)