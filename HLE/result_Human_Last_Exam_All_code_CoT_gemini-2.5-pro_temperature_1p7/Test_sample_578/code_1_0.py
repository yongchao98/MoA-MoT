import math

def calculate_d_n(n):
    """
    Calculates the dimension d_n, which is the (n-1)-th Catalan number.
    The k-th Catalan number is C_k = (1/(k+1)) * C(2k, k).
    So, d_n = C_{n-1} = (1/n) * C(2n-2, n-1).
    """
    if n == 1:
        # C_0 = 1
        return 1
    
    # Use math.comb for combination C(n,k), which is n! / (k! * (n-k)!)
    # math.comb(n, k) is efficient and avoids large intermediate numbers.
    try:
        combinations = math.comb(2 * n - 2, n - 1)
        # The division is guaranteed to be exact.
        return combinations // n
    except ValueError:
        return 0

def solve_product():
    """
    Calculates the product of d_n for n from 1 to 8.
    """
    d_values = []
    product = 1
    
    print("Calculating the dimension d_n for n = 1 to 8:")
    for n in range(1, 9):
        dn = calculate_d_n(n)
        d_values.append(dn)
        product *= dn
        print(f"d_{n} = {dn}")
        
    print("\nThe final product is calculated as follows:")
    
    # Create the full equation string, showing each number
    equation_str = " * ".join(map(str, d_values))
    
    print(f"Product = {equation_str} = {product}")

# Execute the solution
solve_product()