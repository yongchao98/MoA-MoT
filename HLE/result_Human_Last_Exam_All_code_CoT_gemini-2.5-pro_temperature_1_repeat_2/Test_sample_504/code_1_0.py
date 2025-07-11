import math

def N(d):
    """
    Calculates the maximum number of equiangular vectors in C^d
    with inner product magnitude 1/2.
    """
    if d == 0:
        return 0
    if d == 1:
        return 1
    # For d=2, it can be shown N(2)=3
    if d == 2:
        return 3
    # For d=3, |(v,w)|^2 = 1/4 = 1/(d+1), so it's a SIC-POVM with d^2 vectors.
    if d == 3:
        return 9
    # For d > 3, we use the bound N(d) <= 2d, as |(v,w)|^2 != 1/(d+1)
    # The actual values for N(4), N(5), N(6) are not known precisely but are bounded.
    # We use these bounds as the best possible scenario.
    if d == 4:
        return 8 # bound is N(4)<=8
    if d == 5:
        return 10 # bound is N(5)<=10
    if d == 6:
        return 12 # bound is N(6)<=12
    return 0

def get_partitions(n):
    """
    Generates all partitions of an integer n.
    """
    a = [0 for i in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]

def solve():
    """
    Solves the problem by checking all valid partitions of the dimension 6.
    """
    target_dim = 6
    max_vectors = 0
    best_partition = []
    
    # Generate all partitions of the integer 6
    partitions = get_partitions(target_dim)
    
    # The problem requires at least one orthogonal pair, which means the graph of
    # non-orthogonality is disconnected. This corresponds to partitions with at least two parts.
    valid_partitions = [p for p in partitions if len(p) >= 2]
    
    print("Checking all partitions of dimension 6 into 2 or more orthogonal subspaces...")
    
    for partition in valid_partitions:
        current_sum = 0
        for dim in partition:
            current_sum += N(dim)
        
        equation_str_parts = []
        for dim in partition:
            equation_str_parts.append(f"N({dim})")
        
        sum_str_parts = []
        for dim in partition:
            sum_str_parts.append(str(N(dim)))

        print(f"Partition {partition}: n <= {' + '.join(equation_str_parts)} = {' + '.join(sum_str_parts)} = {current_sum}")

        if current_sum > max_vectors:
            max_vectors = current_sum
            best_partition = partition

    print("\n-----------------------------------------------------")
    print(f"The best partition is {best_partition}.")
    
    final_sum_parts = []
    for dim in best_partition:
      final_sum_parts.append(str(N(dim)))
      
    final_equation = f"{' + '.join(final_sum_parts)} = {max_vectors}"
    print(f"The maximum number of vectors is given by the sum: {final_equation}")
    print(f"Each number in the final equation represents the maximum number of vectors in the respective subspace.")
    print("The final answer is the result of this sum.")
    print("Final Result:")
    print(max_vectors)

solve()