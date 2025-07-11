import math

def solve():
    """
    Calculates the product of the dimensions of B_n-invariant subspaces
    of Khovanov homology for T(n,n) torus links, for n=1 to 8.
    """
    
    # p(k) is the number of partitions of integer k. We only need p(0) through p(3).
    # p(0)=1, p(1)=1, p(2)=2, p(3)=3.
    partition_numbers = [1, 1, 2, 3]

    # a_n = d_n / 2. We will find a recurrence for a_n.
    # The sequence a_n starts a_1=1, a_2=2, a_3=2.
    # For n>=4, a_n = a_{n-1} + p(floor((n-2)/2)).
    a = [0] * 9  # Use 1-based indexing
    a[1] = 1
    a[2] = 2
    a[3] = 2
    
    for n in range(4, 9):
        p_index = (n - 2) // 2
        a[n] = a[n - 1] + partition_numbers[p_index]
        
    # d_n = 2 * a_n
    d = [0] * 9 # Use 1-based indexing
    for n in range(1, 9):
        d[n] = 2 * a[n]

    # Calculate the product
    product = 1
    for n in range(1, 9):
        product *= d[n]

    # Format the output string
    equation_parts = [str(d[n]) for n in range(1, 9)]
    equation_str = " * ".join(equation_parts)
    
    print(f"{equation_str} = {product}")
    print(f"The product is: {product}")

solve()