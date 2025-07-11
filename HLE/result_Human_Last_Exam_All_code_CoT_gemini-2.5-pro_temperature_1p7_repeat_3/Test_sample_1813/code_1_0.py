import math

def compute_generalized_markov_cf(p, q):
    """
    Computes the continued fraction of the generalized Markov number m_{p/q}.
    The process is:
    1. Determine the value of m_{p/q}. For r=4/7, this is sqrt(136)/5.
    2. Compute the continued fraction of this value.
    This function implements step 2 for r=4/7.
    """
    # The number is sqrt(136)/5.
    # We use the standard algorithm for quadratic surds (a*sqrt(D)+b)/c.
    # To use the integer version, we convert sqrt(136)/5 to
    # (sqrt(136*25))/25 = sqrt(3400)/25.
    D = 3400
    P = 0
    Q = 25

    # Store (P,Q) pairs to detect the repeating cycle
    history = {}
    coefficients = []
    
    # The algorithm starts
    while (P, Q) not in history:
        history[(P, Q)] = len(coefficients)
        
        # Calculate the next coefficient
        a = math.floor((P + math.sqrt(D)) / Q)
        coefficients.append(a)
        
        # Update P and Q for the next iteration
        P_next = a * Q - P
        Q_next = (D - P_next**2) / Q
        
        P = P_next
        Q = int(Q_next)
    
    # Cycle is found
    start_of_period = history[(P, Q)]
    non_periodic_part = coefficients[:start_of_period]
    periodic_part = coefficients[start_of_period:]

    return non_periodic_part, periodic_part

def main():
    """
    Main function to compute and print the continued fraction for m_{4/7}.
    """
    p, q = 4, 7
    
    # The value of m_{4/7} is sqrt(136)/5.
    D = 136
    val_str = f"sqrt({D})/5"
    
    non_periodic, periodic = compute_generalized_markov_cf(p, q)
    
    # Print the result
    print(f"The continued fraction for the generalized Markov number m_{p}/{q} = {val_str} is:")
    a0 = non_periodic[0]
    period_str = ", ".join(map(str, periodic))
    print(f"[{a0}; ({period_str})]")
    
    print("\nThe equation form is:")
    
    # Construct the equation string
    # We'll show the first few terms including one full period.
    terms_to_show = len(periodic) + 1 
    all_coeffs = [non_periodic[0]] + periodic * 2 # get enough terms
    
    # Print the first coefficient
    print(f"{val_str} = {all_coeffs[0]} + 1 / (", end="")
    
    # Print the rest of the equation nesting parentheses
    for i in range(1, terms_to_show):
        print(f"{all_coeffs[i]} + 1 / (", end="")
    
    print("...)))" + ")" * (terms_to_show - 1))

if __name__ == '__main__':
    main()