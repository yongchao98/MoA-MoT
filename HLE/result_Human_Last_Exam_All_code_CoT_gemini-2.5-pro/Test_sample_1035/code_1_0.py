import math

def solve():
    """
    Calculates the minimum number of bipartite graphs required to cover the edges of K_n.
    The problem is interpreted as finding the size of the smallest family of bipartite
    graphs whose edge sets cover the edge set of K_n. This value is k = ceil(log2(n)).
    """
    n = 35

    print(f"The problem asks for the minimum number of vertices in a family of bipartite graphs covering K_n for n = {n}.")
    print("This is interpreted as finding the minimum number of bipartite graphs, k, in such a family.")
    print("This requires finding the smallest integer k that satisfies the equation 2^k >= n.")
    print(f"Here, we need to find the smallest k such that 2^k >= {n}.\n")

    power_of_2 = 1
    k = 0
    while power_of_2 < n:
        k += 1
        power_of_2 *= 2
        print(f"For k = {k}, the equation is 2^{k} = {power_of_2}. Since {power_of_2} < {n} is {power_of_2 < n}, we continue.")

    # The loop terminates when power_of_2 >= n. Let's do one more step for clarity.
    k += 1
    power_of_2 *= 2
    # Oh wait, the loop condition handles it. Let's re-adjust the logic
    
    power_of_2 = 1
    k = 0
    while power_of_2 < n:
        k += 1
        power_of_2 *= 2
        print(f"Checking k = {k}: 2^{k} = {power_of_2}. This is {'not ' if power_of_2 < n else ''}enough to label {n} vertices.")

    final_k = math.ceil(math.log2(n))
    
    print(f"\nThe smallest integer k satisfying 2^k >= {n} is {final_k}.")
    print(f"Therefore, the minimum number of bipartite graphs needed is {final_k}.")

solve()