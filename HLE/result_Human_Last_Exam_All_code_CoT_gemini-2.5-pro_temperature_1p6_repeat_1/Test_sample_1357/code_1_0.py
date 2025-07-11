import collections

def find_and_count_partitions(g):
    """
    Finds and counts the number of partitions of an integer g.
    This corresponds to the number of stable reduction types for a genus g curve
    with a good reduction Jacobian.
    """
    # We use a recursive helper function to find the partitions
    partitions_list = []
    def find_partitions_recursive(target, max_val, current_partition):
        if target == 0:
            partitions_list.append(list(current_partition))
            return
        
        for i in range(min(max_val, target), 0, -1):
            current_partition.append(i)
            find_partitions_recursive(target - i, i, current_partition)
            current_partition.pop()

    find_partitions_recursive(g, g, [])

    print(f"The problem reduces to finding the number of partitions of the genus g = {g}.")
    print("Each partition represents a distinct type of stable reduction, where the numbers in the partition are the genera of the smooth components.")
    print("-" * 20)
    
    # Print each partition
    print("The possible partitions of the genus are:")
    for p in partitions_list:
        print(f"  {p}")
    
    # Print the final calculation as a sum
    count = len(partitions_list)
    sum_string = " + ".join(["1"] * count)
    print("\nCounting these types gives the total number:")
    print(f"Total = {sum_string} = {count}")
    
# The genus of the curve is 4
genus = 4
find_and_count_partitions(genus)