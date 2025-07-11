import math

def partition_function_memo(n, memo):
    """
    Calculates the partition number p(n) using memoization.
    """
    if n in memo:
        return memo[n]
    if n < 0:
        return 0
    if n == 0:
        return 1
    
    val = 0
    k = 1
    while True:
        pentagonal1 = k * (3 * k - 1) // 2
        pentagonal2 = k * (3 * k + 1) // 2
        
        if pentagonal1 > n and pentagonal2 > n:
            break
        
        term1 = partition_function_memo(n - pentagonal1, memo)
        term2 = partition_function_memo(n - pentagonal2, memo)

        if k % 2 == 1:
            val += term1 + term2
        else:
            val -= term1 + term2
        
        k += 1
        
    memo[n] = val
    return val

def calculate_dn(n, p_memo):
    """
    Calculates d_n based on the derived formula.
    """
    if n == 1:
        # Special case for n=1, T(1,1) is the unknot.
        # d_1 is the dimension of Kh(unknot), which is 2.
        return 2
    
    # For n>=2, d_n = sum_{k=0}^{floor(n/2)} p(n-2k)
    total = 0
    for k in range(math.floor(n / 2) + 1):
        total += partition_function_memo(n - 2 * k, p_memo)
    return total

def main():
    """
    Main function to calculate the product and print the results.
    """
    # Pre-calculate partition numbers up to 8
    p_memo = {}
    for i in range(9):
        partition_function_memo(i, p_memo)
        
    d_values = []
    for n in range(1, 9):
        d_values.append(calculate_dn(n, p_memo))
        
    print("The values of d_n for n=1 to 8 are:")
    for i, d_val in enumerate(d_values):
        print(f"d_{i+1} = {d_val}")
        
    product = 1
    for val in d_values:
        product *= val
        
    # Format the equation string
    equation_str = " * ".join(map(str, d_values))
    
    print("\nThe product is:")
    print(f"{equation_str} = {product}")

if __name__ == "__main__":
    main()
