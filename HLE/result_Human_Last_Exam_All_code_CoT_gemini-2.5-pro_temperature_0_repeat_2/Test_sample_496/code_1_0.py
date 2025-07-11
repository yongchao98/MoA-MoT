def solve_total_rank():
    """
    Calculates the total rank of the SO(4)-equivariant cohomology ring of SO(4) \ X
    in degrees up to 100.

    The problem reduces to calculating the sum of the ranks of the cohomology
    of the classifying space BSO(4). The rational cohomology of BSO(4) is a
    polynomial ring in two generators of degree 4. Its Poincare series is
    P(t) = 1 / (1 - t^4)^2 = sum_{n=0 to inf} (n+1) * t^(4n).

    The rank of the cohomology in degree k=4n is (n+1), and 0 otherwise.
    We need to find the total rank for degrees k <= 100.
    This means we sum the ranks for all n such that 4n <= 100, i.e., n <= 25.
    The total rank is the sum of (n+1) for n from 0 to 25.
    """
    max_degree = 100
    # The degrees are of the form 4n
    max_n = max_degree // 4

    total_rank = 0
    equation_terms = []

    # Sum the ranks for n from 0 to 25
    for n in range(max_n + 1):
        # The rank at degree 4n is n+1
        rank_at_4n = n + 1
        total_rank += rank_at_4n
        equation_terms.append(str(rank_at_4n))

    print("The total rank is the sum of the ranks for degrees 4n, where n ranges from 0 to 25.")
    print("The rank at degree 4n is n+1.")
    print("The calculation is the sum of (n+1) for n=0...25, which is the sum of integers from 1 to 26.")
    
    # To avoid an extremely long output, we show a summary of the equation
    if len(equation_terms) > 10:
        start_terms = " + ".join(equation_terms[:5])
        end_terms = " + ".join(equation_terms[-5:])
        equation_str = f"{start_terms} + ... + {end_terms}"
    else:
        equation_str = " + ".join(equation_terms)

    print(f"\nSum: {equation_str}")
    print(f"Total Rank = {total_rank}")

solve_total_rank()