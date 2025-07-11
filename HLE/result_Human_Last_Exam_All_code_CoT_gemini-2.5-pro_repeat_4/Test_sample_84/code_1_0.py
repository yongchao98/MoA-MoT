import math

def solve_alpha():
    """
    This function calculates the value of alpha based on the problem description.

    The asymptotic growth rate of the minimum polynomial degree d_n is given by
    d_n = Theta(n^alpha).

    The polynomial p_n must satisfy:
    - p_n(i) is in [0,1] for i in {1, 2, ..., n^2}
    - p_n(i) is in [2,3] for i in {n^2+1, ..., n^10}

    This problem on discrete points can be shown to be equivalent to a problem
    of separating two continuous intervals:
    I1 = [1, n^2] and I2 = [n^2+1, n^10].

    The lengths of these intervals are approximately:
    L1 ~ n^exp1
    L2 ~ n^exp2
    The gap between them has length g = 1.
    """

    # From the problem statement, the first set of points extends to n^2.
    exp1 = 2
    # The second set of points extends to n^10.
    exp2 = 10

    print(f"The problem involves two sets of points. The first set is defined up to n^{exp1}.")
    print(f"The second set is defined up to n^{exp2}.")
    
    # A key result from approximation theory states that the degree d_n required to
    # separate two intervals of lengths L1 and L2 with a gap of length g is
    # asymptotically proportional to sqrt(L1 * L2) / g.
    
    # For this problem, L1 is proportional to n^exp1, L2 is proportional to n^exp2, and g=1.
    # So, d_n is proportional to sqrt(n^{exp1} * n^{exp2}).

    print(f"The asymptotic degree d_n follows the relation:")
    print(f"d_n ~ sqrt(n^{exp1} * n^{exp2})")
    print(f"d_n ~ sqrt(n^({exp1} + {exp2}))")
    
    # The exponent alpha is therefore (exp1 + exp2) / 2.
    alpha = (exp1 + exp2) / 2.0
    
    print(f"d_n ~ n^(({exp1} + {exp2}) / 2)")
    print(f"d_n ~ n^{alpha}")
    
    print("\nTherefore, the value of alpha is:")
    print(int(alpha))

solve_alpha()