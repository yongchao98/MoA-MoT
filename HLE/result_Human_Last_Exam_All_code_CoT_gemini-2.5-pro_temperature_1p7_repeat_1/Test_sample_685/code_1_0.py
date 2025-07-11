import sys

def solve():
    """
    Determines if the first player has a winning position with a probability > 50%.
    """
    if len(sys.argv) != 3:
        print("Usage: python your_script.name <n> <m>")
        return

    try:
        n = int(sys.argv[1])
        m = int(sys.argv[2])
        if n <= 0 or m <= 0:
            raise ValueError
    except ValueError:
        print("Please provide positive integers for n and m.")
        return

    # As analyzed in the thinking steps, the first player has a winning probability > 50%
    # if and only if the probability of the initial random board being a P-position
    # is strictly less than 50%.

    # This condition holds for all (n, m) pairs except for n=1, m=1.
    # For n=1, m=1, there are two states: [[0]] (P-position) and [[1]] (N-position).
    # P(P-position) = 0.5, which is not < 0.5.
    # For all other cases, the number of P-positions is a very small fraction of the total
    # 2^(n*m) states, so P(P-position) < 0.5.

    # Therefore, the function f(n, m) returns 1 if (n,m) is not (1,1), and 0 otherwise.
    # This logic can be implemented with a simple conditional check.

    if n == 1 and m == 1:
        result = 0
    else:
        result = 1

    print(f"For n = {n}, m = {m}:")
    print(f"The function f(n, m) returns: {result}")
    
solve()
