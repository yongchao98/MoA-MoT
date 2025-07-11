import math

def max_complex_zeros(N):
    """
    Calculates the maximal possible number of complex zeros for the given matrix B(k).

    Args:
        N (int): The dimension of the square matrix.

    Returns:
        int: The maximal number of complex zeros.
    """
    if not isinstance(N, int) or N <= 0:
        raise ValueError("N must be a positive integer.")
    
    # If N is even, the number is N * 2^(N-1)
    # If N is odd, the number is (N-1) * 2^(N-1)
    # This can be written as (N - (N % 2)) * 2**(N-1)
    
    if N % 2 == 0:
        # N is even
        result = N * (2**(N - 1))
    else:
        # N is odd
        result = (N - 1) * (2**(N - 1))
        
    return result

def main():
    """
    Main function to demonstrate the calculation for N from 1 to 5.
    """
    print("Calculating the maximal possible number of complex zeros for N=1 to 5.")
    for n in range(1, 6):
        num_zeros = max_complex_zeros(n)
        print(f"For N = {n}, the maximal number of complex zeros is {num_zeros}.")

if __name__ == "__main__":
    main()
