import math

def solve_complexity():
    """
    This function analyzes the query complexity for the two specified regimes
    and prints the result in the required (a,b,c) format.
    """

    # The complexity is derived from Theta(min(N*logN, N*L/logN)).
    # The choice depends on comparing L with (logN)^2.
    # If L > (logN)^2, complexity is Theta(N*logN).
    # If L <= (logN)^2, complexity is Theta(N*L/logN).

    # For both regimes analyzed below, the resulting complexity is Theta(N*logN).
    # We now convert this to the (a,b,c) format.
    # Q = N*logN
    # The format is sqrt(N^a * (logN)^b * (loglogN)^c).
    # Q^2 = (N*logN)^2 = N^2 * (logN)^2 * (loglogN)^0
    # By comparing Q^2 with the expression under the square root, we get:
    # a = 2
    # b = 2
    # c = 0
    
    a = 2
    b = 2
    c = 0
    
    # Analysis for Regime 1: N = 2^sqrt(L) implies L = (log2 N)^2.
    # Here, L is equal to (log N)^2, so the complexity is Theta(N*L/logN) = Theta(N * (logN)^2 / logN) = Theta(N*logN).
    # This gives (a,b,c) = (2,2,0).
    result1 = (a, b, c)

    # Analysis for Regime 2: N = 2^((log2 L)^2) implies L = 2^sqrt(log2 N).
    # Here, L grows asymptotically faster than (log N)^2.
    # Thus, the complexity is min(N*logN, N*L/logN) = Theta(N*logN).
    # This also gives (a,b,c) = (2,2,0).
    result2 = (a, b, c)

    # The problem asks to output the numbers in the final equation.
    # We will print the two resulting tuples in the specified format.
    print(f"({result1[0]},{result1[1]},{result1[2]}),({result2[0]},{result2[1]},{result2[2]})")

solve_complexity()