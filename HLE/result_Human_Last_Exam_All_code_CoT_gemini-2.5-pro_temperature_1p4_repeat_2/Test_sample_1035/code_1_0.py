import math

def solve_graph_problem():
    """
    Calculates the minimum number of vertices in a family of bipartite graphs
    covering all the edges of the complete graph K_n for n = 35.
    
    The solution is based on the theorem that the minimum value is n * ceil(log2(n)).
    """
    n = 35

    # Step 1: Calculate k = ceil(log2(n))
    # 2^5 = 32 and 2^6 = 64, so log2(35) is between 5 and 6.
    # The ceiling will be 6.
    k = math.ceil(math.log2(n))

    # Step 2: Calculate the result by multiplying n with k
    result = n * k

    # Step 3: Print the explanation and the final calculation,
    # showing each number in the equation.
    print(f"The complete graph is K_n where n = {n}.")
    print("The minimum total number of vertices in a family of bipartite graphs covering K_n is given by the formula: n * ceil(log2(n)).")
    print(f"First, we calculate the value of ceil(log2({n})).")
    print(f"log2({n}) is approximately {math.log2(n):.4f}")
    print(f"The ceiling of {math.log2(n):.4f} is {k}.")
    print("\nNow, we compute the final result:")
    print(f"{n} * {k} = {result}")

solve_graph_problem()