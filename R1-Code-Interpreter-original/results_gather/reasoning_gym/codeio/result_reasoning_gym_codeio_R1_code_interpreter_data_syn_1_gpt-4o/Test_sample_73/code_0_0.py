def partition_count(n):
    # Create a list to store partition counts
    partitions = [0] * (n + 1)
    partitions[0] = 1  # Base case: there's one way to partition 0

    # Calculate partition numbers using dynamic programming
    for i in range(1, n + 1):
        for j in range(i, n + 1):
            partitions[j] += partitions[j - i]

    return partitions[n]

# Find the integer whose partition count is 89134
target_partition_count = 89134
n = 1
while True:
    if partition_count(n) == target_partition_count:
        print(n)
        break
    n += 1