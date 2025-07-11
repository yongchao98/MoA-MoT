import math
import argparse

def solve():
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.
    The formula is based on a direct construction of an Algebraic Branching Program.
    """
    parser = argparse.ArgumentParser(description="Calculate matrix product complexity.")
    parser.add_argument('--n', type=int, default=5, help='Number of choices for each index.')
    parser.add_argument('--d', type=int, default=3, help='Number of variables in each term (must be odd).')
    args = parser.parse_args()

    n = args.n
    d = args.d

    if not (2 <= d <= n):
        print(f"Error: constraints not met. Required 2 <= d <= n, but got n={n}, d={d}.")
        return

    # The problem states d is odd, but the construction works for any d.
    # We can add a check if we want to be strict.
    if d % 2 == 0:
        print(f"Warning: The problem specified that d ({d}) should be odd.")

    m = []
    for i in range(1, d):
        try:
            # m_i = i * C(n, i)
            m_i = i * math.comb(n, i)
            m.append(m_i)
        except ValueError:
            print(f"Error: n={n}, i={i} is not a valid combination for math.comb (n must be >= k).")
            return

    complexity = 2 + sum(m)
    
    equation_str = " + ".join(map(str, m))
    
    print(f"For n={n} and d={d}:")
    print(f"The dimensions of the intermediate matrices are: {m}")
    print(f"The complexity is {complexity}.")
    print("The equation for the complexity is:")
    print(f"{complexity} = 2 + {equation_str}")


if __name__ == "__main__":
    solve()
