import math

def calculate_ranks():
    """
    Calculates the total rank of the equivariant cohomology ring A up to degree 100.
    The total rank is a sum of two parts, based on the analysis of a long exact sequence.
    """

    # Part 1: Sum of ranks for H_{SO(4)}^*(SO(4))
    # The ring is Q[c1, c2] with deg(c1)=deg(c2)=2.
    # The rank of the k-th degree component (where k=2d) is the number of monomials
    # of degree d in 2 variables, which is d + 1.
    # We sum the ranks for k from 0 to 100, which corresponds to d from 0 to 50.
    sum1 = 0
    for d in range(51): # d from 0 to 50
        sum1 += (d + 1)

    # Part 2: Sum of ranks for H_{SO(4)}^*(X) shifted by -2.
    # The ring is Q[x, p] with deg(x)=2, deg(p)=4.
    # We need the rank of H_G^{k-2}(X), so we are summing ranks for degrees m = k-2.
    # m ranges from -2 to 98. We only consider non-negative even m.
    # The rank for degree m=2d is the number of non-negative integer solutions to i + 2j = d.
    # This number is floor(d/2) + 1.
    # m goes from 0 to 98 (even), so d = m/2 goes from 0 to 49.
    sum2 = 0
    for d_prime in range(50): # d_prime from 0 to 49
        sum2 += (d_prime // 2 + 1)
        
    total_rank = sum1 + sum2

    # The problem asks to output each number in the final equation.
    # The equation is: Total Rank = Part 1 + Part 2
    print("Value for Part 1 (from H_{SO(4)}^*(SO(4))):")
    print(sum1)
    print("Value for Part 2 (from H_{SO(4)}^*(X) term):")
    print(sum2)
    print("The total rank is the sum of these two parts:")
    print(total_rank)

calculate_ranks()