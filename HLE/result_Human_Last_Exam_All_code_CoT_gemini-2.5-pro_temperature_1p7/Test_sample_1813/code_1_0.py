import math

def compute_markov_cf():
    """
    Computes the continued fraction associated with the generalized Markov number m_{4/7}.
    """
    # Step 1: Define the rational number
    p_rat, q_rat = 4, 7

    # Step 2: Compute the continued fraction of p/q
    def get_rational_cf(p, q):
        coeffs = []
        temp_p, temp_q = p, q
        while temp_q != 0:
            a = temp_p // temp_q
            coeffs.append(a)
            temp_p, temp_q = temp_q, temp_p % temp_q
        return coeffs

    cf_rat = get_rational_cf(p_rat, q_rat)

    # Step 3: Construct the matrix K
    def multiply_matrices(A, B):
        C = [[0, 0], [0, 0]]
        C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0]
        C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1]
        C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0]
        C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1]
        return C

    M = lambda a: [[a, 1], [1, 0]]
    K = [[1, 0], [0, 1]]  # Start with identity matrix
    for a in cf_rat:
        K = multiply_matrices(K, M(a))

    # Step 4: Calculate the value of the generalized Markov number m_{4/7}
    trace_K = K[0][0] + K[1][1]
    det_K = K[0][0] * K[1][1] - K[0][1] * K[1][0]
    D_m = trace_K**2 - 4 * det_K

    # Step 5: Compute the continued fraction for m_{4/7} = sqrt(D_m)/q_rat
    # To use the integer algorithm, we express m as (sqrt(D_m*q_rat^2) + 0) / q_rat^2
    D = D_m * q_rat**2
    P = 0
    Q = q_rat**2

    coeffs_m = []
    history = {}  # Store (P, Q) states to detect the cycle

    # Loop to find the coefficients and the period
    for i in range(50):  # Safety limit
        state = (P, Q)
        if state in history:
            pre_period_start_index = history[state]
            pre_period = coeffs_m[:pre_period_start_index]
            period = coeffs_m[pre_period_start_index:]
            break
        history[state] = i
        
        a = int((math.sqrt(D) + P) / Q)
        coeffs_m.append(a)
        
        P_next = a * Q - P
        Q_next = (D - P_next**2) // Q
        
        P, Q = P_next, Q_next
    else: # Should not be reached for quadratic irrationals
        pre_period = coeffs_m
        period = []

    # Step 6: Format and print the final result in a human-readable equation.
    print(f"The continued fraction associated with the generalized Markov number m_{p_rat}/{q_rat} is:")
    
    # Construct the output string
    output_str = f"m_{p_rat}/{q_rat} = [{pre_period[0]};"

    # Add non-repeating part (if any)
    if len(pre_period) > 1:
        output_str += " " + ", ".join(map(str, pre_period[1:]))
        if period:
            output_str += ","
    
    # Add repeating part using overline notation
    if period:
        output_str += " " + "overline(" + ", ".join(map(str, period)) + ")"

    output_str += "]"
    print(output_str)

if __name__ == '__main__':
    compute_markov_cf()