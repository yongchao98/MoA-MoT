def find_integer_partitions(n):
    """
    Generates all unique partitions of an integer n.
    Each partition is a list of integers that sum to n.
    The function uses a recursive helper to find the partitions.
    """
    partitions_list = []
    
    # A recursive helper function to find partitions
    # target: the remaining sum we need to achieve
    # max_val: the maximum value a part can have (to ensure non-increasing order)
    # current_partition: the list of parts found so far
    def generate(target, max_val, current_partition):
        if target == 0:
            # When the target is 0, a valid partition has been found
            partitions_list.append(list(current_partition))
            return

        # Iterate from min(target, max_val) down to 1
        # This ensures parts are in non-increasing order, avoiding duplicates like [3,1] and [1,3]
        for i in range(min(target, max_val), 0, -1):
            current_partition.append(i)
            generate(target - i, i, current_partition)
            # Backtrack to explore other possibilities
            current_partition.pop()

    generate(n, n, [])
    return partitions_list

def solve_problem():
    """
    Main function to solve the user's question.
    """
    genus = 4
    
    print(f"The number of types of stable reductions of a genus {genus} curve")
    print("with good reduction Jacobian corresponds to the number of integer partitions of the genus.\n")
    
    partitions = find_integer_partitions(genus)
    
    print(f"The {len(partitions)} distinct partitions of {genus} are:")
    for p in partitions:
        print(" -> " + " + ".join(map(str, p)))

    count = len(partitions)
    
    # Build and print the final equation as requested
    equation_parts = ['1'] * count
    final_equation = " + ".join(equation_parts)
    
    print(f"\nEach partition represents one type. Summing these types gives the total:")
    print(f"{final_equation} = {count}")
    
    # Final answer in the required format
    print(f"\n<<<{count}>>>")

# Run the solver
solve_problem()