import math

def calculate_limit_ratio():
    """
    Calculates the ratio of the expected number of remaining numbers to n
    and demonstrates its convergence to the analytical limit.
    """
    # E_n: Expected number of remaining items for a sequence of length n.
    # S_n: Sum of E_k from k=0 to n.
    # We start with E_0 = 0 and E_1 = 1.
    E = [0.0, 1.0]
    S = [0.0, 1.0]  # S[n] stores sum(E_k for k in 0..n)

    max_n = 10000

    print("Calculating the ratio E_n/n for increasing values of n:")
    
    # We calculate E_n using the recurrence E_n = (2 / (n-1)) * S_{n-2}
    for n in range(2, max_n + 1):
        # S_{n-2} is stored at index n-2 in the list S
        e_n = (2.0 / (n - 1)) * S[n - 2]
        E.append(e_n)
        
        # Update the cumulative sum S_n = S_{n-1} + E_n
        s_n = S[n - 1] + e_n
        S.append(s_n)

        # Print the ratio at specific intervals to observe convergence
        if n in [10, 100, 1000, 10000]:
            ratio = E[n] / n
            print(f"For n = {n:5d}, the ratio E_n/n is: {ratio:.8f}")

    # The analytical limit is 1/4 - 3/(4*e^2)
    limit_value = 1.0/4.0 - 3.0/(4.0 * math.exp(2))
    
    print("\n--------------------------------------------------")
    print("The exact limit as n approaches infinity is given by the equation:")
    print("L = 1/4 - 3/(4 * e^2)")
    print("--------------------------------------------------\n")
    
    print("The numbers and symbols in the final equation are:")
    print("1")
    print("4")
    print("3")
    print("4")
    print("e (Euler's number)")
    print("2")
    
    print(f"\nThe numerical value of this limit is approximately: {limit_value:.8f}")

calculate_limit_ratio()