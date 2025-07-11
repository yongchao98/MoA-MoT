def solve():
    """
    Calculates the total rank of the SO(4)-equivariant cohomology ring of SO(4)\\X 
    for degrees up to 100, based on its Poincaré series.

    The Poincaré series for the ranks of the cohomology ring A is:
    P_A(t) = (t^3 + t^4 + t^6 + t^7) / (1 - t^4)^2
    
    The rank in a given degree k, d_k, is the coefficient of t^k in the series expansion.
    The total rank up to degree 100 is the sum of d_k for k from 0 to 100.
    """
    
    limit = 100
    
    # The offsets correspond to the powers in the numerator of the Poincaré series.
    offsets = [3, 4, 6, 7]
    
    results = []
    
    # For each term t^j in the numerator, we calculate the sum of the corresponding
    # coefficients in the final series expansion.
    for offset in offsets:
        # We need to sum coefficients for terms of the form t^(4n+offset) <= t^limit.
        # This means 4n <= limit - offset, so n runs up to floor((limit - offset)/4).
        # The coefficient for t^(4n) in 1/(1-t^4)^2 is (n+1).
        # We need to sum (n+1) for n from 0 to n_max.
        if limit >= offset:
            n_max = (limit - offset) // 4
            # The sum of integers from 1 to (n_max + 1) is (n_max + 1) * (n_max + 2) / 2.
            current_sum = (n_max + 1) * (n_max + 2) // 2
            results.append(current_sum)
        else:
            # If offset > limit, this term doesn't contribute to the sum.
            results.append(0)

    total_rank = sum(results)
    
    # As requested, output the numbers in the final equation.
    # These numbers are the total rank contributions from each term in the numerator.
    print(f"{results[0]} + {results[1]} + {results[2]} + {results[3]} = {total_rank}")

solve()