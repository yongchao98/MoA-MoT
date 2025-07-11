import math

def main():
    """
    Calculates the product of the dimensions of the B_n-invariant
    subspaces of Kh(T(n,n); Q) for n=1 to 8.
    """
    d_values = []
    # Calculate d_n for n from 1 to 8
    for n in range(1, 9):
        # The formula for d_n is 2 * C(2n-2, n-1)
        # where C is the binomial coefficient.
        try:
            dn = 2 * math.comb(2 * n - 2, n - 1)
            d_values.append(dn)
        except ValueError:
            # math.comb throws ValueError for n < k
            # This case shouldn't be reached for n>=1
            pass

    # Calculate the product of all d_n values
    product = math.prod(d_values)

    # Prepare the output string for the equation
    equation_str = " * ".join(map(str, d_values))
    
    # Print the results as requested
    print(equation_str + " = " + str(product))

if __name__ == "__main__":
    main()