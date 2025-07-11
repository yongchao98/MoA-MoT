def solve():
    """
    This function analyzes the graph peeling process and determines the asymptotic bounds
    on the maximum number of steps for three different classes of forests.

    The analysis reveals that a hierarchical graph structure can be built to make the
    process last for a number of steps K related to the number of vertices n and the
    maximum degree Delta by the formula: K is approximately proportional to log(n)/log(Delta).

    To find the maximum number of steps for a class of graphs, we need to maximize this
    quantity over all graphs in the class. The number of steps K is maximized when
    Delta is as small as possible. A small constant Delta (e.g., Delta=3) can be
    achieved for this construction.

    1. For forests with max degree at most sqrt(log n), we can choose a graph with a
       small constant degree (which is <= sqrt(log n)). This gives a maximum number
       of steps K = Theta(log n). This corresponds to option 8.

    2. For forests with max degree at most log n, we can similarly choose a graph with
       a small constant degree. This gives a maximum number of steps K = Theta(log n).
       This corresponds to option 8.

    3. For any forest, we are free to choose the graph. To maximize K, we again pick
       a small constant degree. This gives K = Theta(log n). This corresponds to option 8.
    """

    # The i-th digit corresponds to f_i(n)
    # 1: O(1)
    # 2: omega(1) but O(log^0.1(log n))
    # 3: omega(log^0.1(log n)) but O(log^0.9(log n))
    # 4: omega(log^0.9(log n)) but O(log(log n))
    # 5: omega(log(log n)) but 2^O(sqrt(log(log n)))
    # 6: 2^omega(sqrt(log(log n))) but O(log^0.9 n)
    # 7: omega(log^0.9 n) but o(log n)
    # 8: Theta(log n)
    # 9: omega(log n)

    f1_category = 8
    f2_category = 8
    f3_category = 8

    result = f"{f1_category}{f2_category}{f3_category}"
    print(result)

solve()
# The final answer is the three-digit number.
# The thinking process leads to 888.
# Printing the final answer in the desired format.
# Let me format my thinking into the final response and put the code block at the end.
final_answer = 888
print(f'<<<{final_answer}>>>')
