import sys

def calculate_sum_approximation(n):
    """
    Calculates the approximation for the sum S(n) for a given positive integer n
    using the formula: n^2/2 + 1/120 + 1/(252*n).

    This function prints the formula with the value of n substituted,
    and the final numerical result.
    """
    # The approximation formula is derived from the Euler-Maclaurin expansion.
    # The coefficients/divisors in the formula are:
    term1_den = 2
    
    term2_num = 1
    term2_den = 120

    term3_num = 1
    term3_den = 252
    
    # Calculate each term of the formula
    val1 = (n**2) / term1_den
    val2 = term2_num / term2_den
    val3 = term3_num / (term3_den * n)
    
    result = val1 + val2 + val3

    # Outputting the equation with each number explicitly shown as requested
    print(f"The approximation for n={n} is calculated as:")
    print(f"({n}^2 / {term1_den}) + ({term2_num} / {term2_den}) + ({term3_num} / ({term3_den} * {n}))")
    print(f"= {val1} + {val2} + {val3}")
    print(f"= {result}")

def main():
    """
    Main function to parse command-line arguments and run the calculation.
    """
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} <n>")
        print("Where <n> is a positive integer.")
        sys.exit(1)

    try:
        n_str = sys.argv[1]
        n = int(n_str)
        if n <= 0:
            raise ValueError("Input must be a positive integer.")
    except ValueError:
        print("Error: Please provide a positive integer for n.")
        sys.exit(1)
        
    calculate_sum_approximation(n)

if __name__ == "__main__":
    main()