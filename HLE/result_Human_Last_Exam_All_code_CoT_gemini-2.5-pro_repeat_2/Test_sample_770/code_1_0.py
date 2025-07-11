def get_partitions(n):
    """
    Generate all unique partitions of an integer n in lexicographical order.
    Example: get_partitions(4) -> [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
    """
    # Memoization cache to store results for subproblems
    cache = {}
    def partitions_recursive(target, max_val):
        # If result is cached, return it
        if (target, max_val) in cache:
            return cache[(target, max_val)]
        
        # Base case: a partition of 0 is an empty list
        if target == 0:
            return [[]]
        
        # Invalid state
        if target < 0 or max_val == 0:
            return []
        
        # Combine partitions that use max_val and partitions that don't
        
        # Case 1: Use at least one part of size max_val
        partitions_with_max_val = []
        if target >= max_val:
            sub_partitions = partitions_recursive(target - max_val, max_val)
            for p in sub_partitions:
                partitions_with_max_val.append([max_val] + p)
        
        # Case 2: Do not use any part of size max_val
        partitions_without_max_val = partitions_recursive(target, max_val - 1)
        
        result = partitions_with_max_val + partitions_without_max_val
        cache[(target, max_val)] = result
        return result

    return partitions_recursive(n, n)

def solve_problem():
    """
    Calculates the rank of H^2_c(Y, Q) based on the McKay correspondence.
    """
    n = 5  # For the alternating group A_5

    print("The problem is to find the rank of H^2_c(Y, Q) for a crepant resolution Y of X = C^3/I,")
    print("where I is the icosahedral group, isomorphic to A_5.")
    print("\nAccording to the 3D McKay Correspondence, this rank is equal to the number of non-trivial")
    print("conjugacy classes of the group I.")
    print("\nWe will calculate this by finding the total number of conjugacy classes of A_5 and subtracting 1.")

    print(f"\n--- Calculating Conjugacy Classes of A_{n} ---")
    print("A conjugacy class of S_n, represented by a partition p, consists of even permutations if n - len(p) is even.")
    print("An even class splits into two classes in A_n if and only if all parts of p are odd and distinct.")

    all_partitions = get_partitions(n)
    
    total_classes_A_n = 0
    class_contributions = []

    print(f"\nAnalyzing partitions of {n}:")
    for p in all_partitions:
        num_cycles = len(p)
        # A permutation's sign is (-1)^(n - num_cycles)
        is_even = (n - num_cycles) % 2 == 0
        
        print(f"\nPartition p = {p}:")
        if not is_even:
            print(f"  - Corresponds to odd permutations (n - #cycles = {n} - {num_cycles} = {n - num_cycles}).")
            continue
        
        print(f"  - Corresponds to even permutations (n - #cycles = {n} - {num_cycles} = {n - num_cycles}).")
        
        # The identity class is special but follows the general rule.
        if p == [1] * n:
             print("  - This is the identity class. It does not split.")
             contribution = 1
        else:
            # Check the splitting condition
            all_odd = all(part % 2 != 0 for part in p)
            all_distinct = len(set(p)) == len(p)
            
            if all_odd and all_distinct:
                print("  - All parts are odd and distinct. The class splits.")
                contribution = 2
            else:
                print("  - The class does not split.")
                contribution = 1
        
        total_classes_A_n += contribution
        if p != [1] * n:
            class_contributions.append(contribution)

    # Add the identity class contribution at the end for clarity in the sum
    class_contributions.append(1)

    print("\n--- Summary ---")
    equation_str = " + ".join(map(str, class_contributions))
    print(f"Total number of conjugacy classes in A_{n} is the sum of contributions from each even partition:")
    print(f"Total classes = {equation_str} = {total_classes_A_n}")
    
    final_answer = total_classes_A_n - 1
    print("\nThe number of non-trivial classes is the total number of classes minus the identity class.")
    print(f"Final Equation: {total_classes_A_n} - 1 = {final_answer}")
    print(f"\nThus, the rank of H^2_c(Y, Q) is {final_answer}.")
    
    return final_answer

if __name__ == "__main__":
    rank = solve_problem()
    # The final answer is returned in the required format.
    print(f"\n<<<{rank}>>>")
