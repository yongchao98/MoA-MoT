import math

def calculate_simplices():
    """
    Calculates the number of n-simplices for N=200, k=13, and n from 0 to 5.
    
    The number of n-simplices in the nerve of the over-category N_.(Z_N)_{k/}
    is given by the number of non-decreasing sequences of length n+1,
    m_0, m_1, ..., m_n, such that k <= m_i <= N for all i.
    This is a combination with repetition problem, and the formula is:
    C( (N - k + 1) + (n + 1) - 1, n + 1) = C(N - k + n + 1, n + 1)
    """
    N = 200
    k = 13
    
    print(f"Calculating the number of n-simplices for N = {N} and k = {k}.")
    print("The formula is: C(N - k + n + 1, n + 1)")
    print("-" * 50)
    
    results = []
    for n in range(6):
        # Parameters for the binomial coefficient C(m, r)
        m = N - k + n + 1
        r = n + 1
        
        # Calculate the number of simplices using math.comb
        num_simplices = math.comb(m, r)
        results.append(num_simplices)
        
        # Print the detailed calculation for each n
        print(f"For n = {n}:")
        print(f"Number of simplices = C({N} - {k} + {n} + 1, {n} + 1) = C({m}, {r}) = {num_simplices}")
        print()
        
    # The final answer is a list of these numbers
    final_answer = ", ".join(map(str, results))
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    calculate_simplices()