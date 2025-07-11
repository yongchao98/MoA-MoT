import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    i = 3
    while i * i <= n:
        if n % i == 0:
            return False
        i += 2
    return True

def get_integer_partitions(n):
    """Generates all unique integer partitions of n."""
    memo = {}
    def find_partitions(k):
        if k in memo:
            return memo[k]
        if k == 0:
            return [[]]
        if k < 0:
            return []
        
        result_set = set()
        for i in range(1, k + 1):
            for sub_partition in find_partitions(k - i):
                new_partition = tuple(sorted([i] + sub_partition))
                result_set.add(new_partition)
        
        # Sort for canonical representation
        result_list = sorted([list(p) for p in result_set])
        memo[k] = result_list
        return result_list
        
    return find_partitions(n)


def describe_filled_groups(q, m):
    """
    Describes the nonabelian filled groups of order 2*q^m.
    
    Args:
        q (int): An odd prime number.
        m (int): A natural number.
    """
    # 1. Theoretical introduction
    print("### Classification of Nonabelian Filled Groups of Order 2*q^m ###\n")
    print(f"For a group G of order 2*q^m (q is an odd prime, m is a natural number), the following conditions determine if it is a nonabelian filled group:")
    print("1. G must be a semidirect product A \u22ca Z_2, where A is the Sylow q-subgroup of order q^m.")
    print("2. The subgroup A must be abelian.")
    print("3. The action of the non-identity element of Z_2 on A must be the inversion map: a \u21a6 a\u207b\u00b9.")
    print("These groups are known as generalized dihedral groups, Dih(A).\n")
    print("The structure of the abelian group A depends on the partitions of m.\n")

    # 2. Find partitions for m
    partitions = get_integer_partitions(m)

    # 3. Print the results for each group
    print(f"### Results for q={q}, m={m} ###\n")
    print(f"The order is 2 * {q}^{m} = {2 * (q**m)}.")
    print(f"There are {len(partitions)} non-isomorphic abelian groups A of order {q**m}, leading to {len(partitions)} distinct filled groups.\n")
    
    for i, p in enumerate(partitions, 1):
        # Build string representation of A
        a_parts_symbolic = [f"Z_(q^{exp})" for exp in p]
        a_parts_numeric = [f"Z_({q**exp})" for exp in p]
        a_str_symbolic = " \u00d7 ".join(a_parts_symbolic).replace("^1", "")
        a_str_numeric = " \u00d7 ".join(a_parts_numeric)

        print(f"--- Group {i} (from partition {p}) ---")
        print(f"  - The abelian group is A = {a_str_symbolic}")
        print(f"    In numbers: A \u2245 {a_str_numeric}")
        
        # Build the presentation of G = Dih(A)
        generators_a = [chr(ord('a') + j) for j in range(len(p))]
        relations = [f"x^2 = 1"]
        
        # Relations for A
        for j, exp in enumerate(p):
            relations.append(f"{generators_a[j]}^{q**exp} = 1")
        
        # Commutativity relations for A
        for j1 in range(len(p)):
            for j2 in range(j1 + 1, len(p)):
                relations.append(f"{generators_a[j1]}{generators_a[j2]} = {generators_a[j2]}{generators_a[j1]}")

        # Inversion relations for the action of x
        for gen in generators_a:
            relations.append(f"x{gen}x = {gen}^-1") # Using x*g*x since x is its own inverse

        presentation = f"G = <{', '.join(generators_a)}, x | {', '.join(relations)}>"
        
        print(f"  - The corresponding filled group is G = Dih(A) with presentation:")
        print(f"    {presentation}\n")


if __name__ == '__main__':
    try:
        q_str = input("Enter an odd prime number q: ")
        q = int(q_str)
        if not is_prime(q) or q == 2:
            print("Error: q must be an odd prime number.")
        else:
            m_str = input("Enter a natural number m (m > 0): ")
            m = int(m_str)
            if m < 1:
                print("Error: m must be a natural number (positive integer).")
            else:
                describe_filled_groups(q, m)
    except ValueError:
        print("Invalid input. Please enter integers.")
    except Exception as e:
        print(f"An error occurred: {e}")
