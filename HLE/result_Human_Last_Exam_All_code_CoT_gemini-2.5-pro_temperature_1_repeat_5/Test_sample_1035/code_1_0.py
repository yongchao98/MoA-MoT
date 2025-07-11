import math

def solve():
    """
    Calculates the minimum number of vertices in a family of bipartite graphs
    covering all the edges of the complete graph K_n.

    This is interpreted as the minimum sum of the sizes of the vertex sets
    of the bipartite graphs in the covering family, a quantity denoted bs(K_n).

    The formula is: bs(K_n) = n*k - max_weight_sum, where:
    k = ceil(log2(n))
    m = 2^k - n
    max_weight_sum is the maximum sum of Hamming weights of m distinct binary
    vectors of length k.
    """
    n = 35

    # Step 1: Calculate k, the minimum number of bipartite graphs.
    k = math.ceil(math.log2(n))

    # Step 2: Calculate m, the number of vectors to select.
    two_k = 2**k
    m = two_k - n

    # Step 3: Calculate the maximum sum of Hamming weights for m vectors of length k.
    # We greedily choose vectors with the most 1s.
    max_weight_sum = 0
    m_remaining = m

    print(f"For n = {n}:")
    print(f"The minimum number of bipartite graphs required is k = ceil(log2({n})) = {k}.")
    print(f"We have {two_k} possible binary vectors of length {k}.")
    print(f"m = {two_k} - {n} = {m}. We need to select {m} vectors with the maximum total Hamming weight.")
    print("\nSelecting vectors greedily based on Hamming weight (w):")

    for w in range(k, -1, -1):
        if m_remaining == 0:
            break
        
        # Number of vectors of length k with weight w
        count_at_weight_w = math.comb(k, w)
        
        # Number of vectors to take from this weight group
        num_to_take = min(m_remaining, count_at_weight_w)
        
        if num_to_take > 0:
            print(f"- Weight w={w}: There are {count_at_weight_w} vectors. We take {num_to_take}.")
            max_weight_sum += num_to_take * w
            m_remaining -= num_to_take
    
    # Step 4: Calculate the final result using the formula.
    result = n * k - max_weight_sum
    
    print(f"\nThe maximum sum of weights for {m} vectors is {max_weight_sum}.")
    print("\nThe final calculation is:")
    print(f"{n} * {k} - {max_weight_sum} = {result}")

solve()
<<<93>>>