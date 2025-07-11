import collections

def solve():
    """
    Solves the problem by finding the optimal partition of the 6-dimensional space.
    """
    
    # N(d) is the maximum number of vectors in C^d with pairwise angle pi/3.
    # This corresponds to |(v,w)| = 1/2.
    # Such sets only exist for d < 4.
    N = {
        1: 1,
        2: 3,
        3: 9
    }
    
    allowed_dims = sorted(list(N.keys()), reverse=True)
    
    memo = {}

    def find_partitions(target_dim):
        """
        Finds all partitions of target_dim into allowed part sizes.
        Returns a list of partitions, where each partition is a list of integers.
        """
        if target_dim in memo:
            return memo[target_dim]
        
        if target_dim == 0:
            return [[]]
        if target_dim < 0:
            return []

        all_partitions = []
        for part in allowed_dims:
            if target_dim >= part:
                sub_partitions = find_partitions(target_dim - part)
                for sub_partition in sub_partitions:
                    new_partition = sorted([part] + sub_partition)
                    if new_partition not in all_partitions:
                        all_partitions.append(new_partition)
        
        memo[target_dim] = all_partitions
        return all_partitions

    dimension = 6
    partitions = find_partitions(dimension)
    
    max_vectors = 0
    best_partition = None
    
    print("Analyzing partitions of dimension 6...")
    print("N(d) is the max number of vectors in C^d with angle pi/3 between any pair.")
    print(f"Known values: N(1)={N[1]}, N(2)={N[2]}, N(3)={N[3]}.\n")

    for p in partitions:
        # The problem states there is at least one pair of orthogonal vectors.
        # A set of equiangular lines cannot contain an orthogonal pair.
        # This implies the set of vectors must span at least two orthogonal subspaces.
        # So, the partition must have at least two parts (m >= 2).
        if len(p) < 2:
            continue
            
        current_vectors = sum(N[dim] for dim in p)
        
        equation_parts = [f"N({dim})" for dim in p]
        sum_parts = [str(N[dim]) for dim in p]
        
        print(f"Partition: {p}")
        print(f"  Total vectors = {' + '.join(equation_parts)} = {' + '.join(sum_parts)} = {current_vectors}")

        if current_vectors > max_vectors:
            max_vectors = current_vectors
            best_partition = p

    print("\n----------------------------------------------------")
    print(f"The best partition is {best_partition}, which splits the space into orthogonal subspaces of dimensions {best_partition}.")
    final_equation_parts = [f"N({dim})" for dim in best_partition]
    final_sum_parts = [str(N[dim]) for dim in best_partition]
    
    print("The final calculation is:")
    print(f"{' + '.join(final_equation_parts)} = {' + '.join(final_sum_parts)} = {max_vectors}")
    print(f"The largest number of vectors is {max_vectors}.")

solve()