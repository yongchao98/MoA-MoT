import math

def solve_ratio_limit(n_max=2000):
    """
    Calculates the expected number of remaining items En and the ratio En/n.
    
    The function computes En based on the recurrence relation:
    (n-1)*En = (n-2)*E(n-1) + 2*E(n-2)
    with base cases E0=0, E1=1, E2=0.
    """
    if n_max < 3:
        raise ValueError("n_max must be at least 3")

    # E stores the expected values E_n
    E = [0.0] * (n_max + 1)
    E[0] = 0
    E[1] = 1
    E[2] = 0

    # Calculate E_n for n from 3 to n_max using the recurrence
    for n in range(3, n_max + 1):
        E[n] = ((n - 2) * E[n - 1] + 2 * E[n - 2]) / (n - 1)

    # The ratio E_n / n for the largest n computed
    ratio = E[n_max] / n_max
    
    # The theoretical limit
    limit = math.exp(-2)

    print(f"The ratio E_n/n for n = {n_max} is: {ratio}")
    print(f"The theoretical limit is e^(-2), which is approximately: {limit}")
    
    # The final equation as requested
    print("\nFinal equation:")
    # The format "output each number in the final equation" is interpreted
    # as showing the components of the limit equation.
    equation_char_by_char = ['e', '^', '(', '-', '2', ')']
    print(' '.join(equation_char_by_char))
    print(f"= {limit}")


solve_ratio_limit()
