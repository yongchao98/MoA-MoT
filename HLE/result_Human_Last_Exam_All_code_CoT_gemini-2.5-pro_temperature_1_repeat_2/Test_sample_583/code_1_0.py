import math

def f(n):
    """
    Calculates f(n) based on the formula derived from the IMO Shortlist 2010 N6 problem.
    Let n-1 = sum(c_i * 2^i) in binary.
    f(n) = 1 + sum(c_i * (i * 2^i + 1))
         = 1 + sum(c_i) + sum(c_i * i * 2^i)
    """
    if n <= 0:
        raise ValueError("n must be a positive integer")
    if n == 1:
        return 1

    m = n - 1
    
    # s2_m corresponds to sum(c_i)
    s2_m = 0
    # sum_term corresponds to sum(c_i * i * 2^i)
    sum_term = 0
    
    i = 0
    temp_m = m
    while temp_m > 0:
        # Check if the i-th bit (c_i) is 1
        if temp_m & 1:
            s2_m += 1
            sum_term += i * (1 << i)  # i * 2^i
        
        temp_m >>= 1
        i += 1
        
    return 1 + s2_m + sum_term

def main():
    """
    Calculates and prints the ratio f(n)/(n*log2(n)) for various large n
    to numerically verify that the limit approaches 1.
    """
    print("This script calculates the value of the expression f(n) / (n * log2(n)) for increasing values of n.")
    print("The theoretical limit as n approaches infinity is 1.")
    print("\nHere are the values for n = 2^k for k from 4 to 24:")
    print("-" * 50)
    print(f"{'n':<12}{'f(n)':<18}{'Ratio f(n)/(n*log2(n))'}")
    print("-" * 50)

    for k in range(4, 25):
        n = 1 << k  # n = 2^k
        
        # Calculate f(n) using the function
        val_f = f(n)
        
        # Calculate the denominator n * log2(n)
        # For n = 2^k, log2(n) is simply k
        denominator = n * k
        
        # Calculate the ratio
        ratio = val_f / denominator
        
        # Print the results in a formatted table
        # Each number in the "final equation" (the ratio calculation) is displayed
        print(f"{n:<12d}{val_f:<18d}{ratio:.10f}")
        
    print("-" * 50)

if __name__ == "__main__":
    main()