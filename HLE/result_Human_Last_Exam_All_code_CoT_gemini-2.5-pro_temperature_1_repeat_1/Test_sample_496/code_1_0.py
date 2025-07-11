import sys

def solve():
    """
    This function calculates the total rank of the equivariant cohomology ring A
    for degrees up to 100.

    The Poincare series for the ranks of A is P(t) = sum_{n=0 to inf} (n+1) * t^(4n+2).
    We need to sum the ranks for all degrees k <= 100.
    The degrees are k = 4n+2.
    So, we need to find the maximum n such that 4n+2 <= 100.
    4n <= 98  =>  n <= 24.5.
    So, n ranges from 0 to 24.
    The total rank is the sum of (n+1) for n from 0 to 24, which is the sum of integers from 1 to 25.
    """

    # The maximum value for n is 24. The sum is of (n+1).
    limit = 25
    
    # Generate the list of numbers to be summed
    numbers_to_sum = [str(i) for i in range(1, limit + 1)]
    
    # Format the equation string
    equation_str = " + ".join(numbers_to_sum)
    
    # Print the full equation
    print(f"Total Rank = {equation_str}")
    
    # Calculate the sum
    total_rank = sum(range(1, limit + 1))
    
    # Print the final result
    print(f"= {total_rank}")

solve()