import math

def get_continued_fraction_for_quadratic_irrational(D, P0, Q0):
    """
    Computes the continued fraction of a quadratic irrational number of the form (P0 + sqrt(D)) / Q0.
    
    The algorithm is based on the standard procedure for such numbers.
    We need to ensure that Q0 divides (D - P0^2).
    
    For alpha = (-13 + sqrt(281))/14, we have P0 = -13, D = 281, Q0 = 14.
    Check: (D - P0^2) / Q0 = (281 - (-13)^2) / 14 = (281 - 169) / 14 = 112 / 14 = 8.
    The condition holds.
    """

    # Check if the number is already negative. If so, P0 has to be handled carefully.
    # Our alpha is positive, so no special handling needed.
    
    P = P0
    Q = Q0
    
    # Check if a0 is negative. 
    # alpha = (-13 + sqrt(281))/14 approx (-13 + 16.76)/14 = 3.76/14 approx 0.26
    # so a0 is 0.
    
    a0 = math.floor((P + math.sqrt(D)) / Q)
    
    coefficients = [a0]
    
    # Store (P, Q) pairs to detect the cycle
    seen_states = {(P, Q)}
    
    while True:
        # Update P and Q for the next step of the fraction
        P_next = coefficients[-1] * Q - P
        Q_next = (D - P_next**2) // Q
        
        P = P_next
        Q = Q_next
        
        # Check if we have entered a cycle
        if (P, Q) in seen_states:
            # Cycle detected. Find the start of the periodic part.
            
            # Recompute from the start to find where the cycle begins
            temp_P, temp_Q = P0, Q0
            non_periodic_part = []
            
            # Get the non-periodic part first
            a_i = math.floor((temp_P + math.sqrt(D)) / temp_Q)
            non_periodic_part.append(a_i)
            
            while (temp_P, temp_Q) != (P, Q):
                 temp_P_next = non_periodic_part[-1] * temp_Q - temp_P
                 temp_Q_next = (D - temp_P_next**2) // temp_Q
                 temp_P, temp_Q = temp_P_next, temp_Q_next
                 a_i = math.floor((temp_P + math.sqrt(D)) / temp_Q)
                 non_periodic_part.append(a_i)
            
            periodic_part = non_periodic_part.pop() # Last element is start of period
            
            return non_periodic_part, [periodic_part] + coefficients[len(non_periodic_part):]

        seen_states.add((P, Q))
        
        # Compute the next coefficient
        a_i = math.floor((P + math.sqrt(D)) / Q)
        coefficients.append(a_i)
        
        # Safety break for very long fractions, though for quadratic irrationals it should always terminate
        if len(coefficients) > 100:
             return "Error: Fraction too long, possible issue.", []


# Parameters for alpha = (-13 + sqrt(281)) / 14
D_val = 281
P0_val = -13
Q0_val = 14

non_periodic, periodic = get_continued_fraction_for_quadratic_irrational(D_val, P0_val, Q0_val)

# Format the output string
periodic_str = ", ".join(map(str, periodic))
print(f"The continued fraction is [{non_periodic[0]}; ({periodic_str})]")
# For clarity, printing each component
print(f"The non-periodic part is: {non_periodic}")
print(f"The periodic part is: {periodic}")
# Output the final requested equation
print(f"m_{4/7} ~ [ {non_periodic[0]}; {periodic[0]}, {periodic[1]}, {periodic[2]}, {periodic[3]}, {periodic[4]}, {periodic[5]}, {periodic[6]}, {periodic[7]}, {periodic[8]}, {periodic[9]}, {periodic[10]}, {periodic[11]}, {periodic[12]}, {periodic[13]}, {periodic[14]}... ]")
