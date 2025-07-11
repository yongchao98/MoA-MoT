import math

def calculate_phi(n):
    """
    Calculates the value of phi(n) based on the derived formula.

    The formula for phi(n) is exp(Tr(Y)), where Tr(Y) has been derived as:
    Tr(Y) = (2 * (n+1) * (n-1)**2) / n
    """
    if not isinstance(n, int) or n < 5:
        print("Error: n must be an integer greater than or equal to 5.")
        return

    # Components of the trace formula's numerator
    c1 = 2
    term1 = n + 1
    term2 = (n - 1)**2
    
    # Numerator and denominator of the trace
    numerator_trace = c1 * term1 * term2
    denominator_trace = n
    
    # Calculate the trace
    trace_val = numerator_trace / denominator_trace
    
    # Calculate phi(n) which is the exponential of the trace
    phi_n = math.exp(trace_val)
    
    # Print the breakdown of the calculation as requested
    print(f"The symbolic formula is: phi(n) = exp( 2 * (n+1) * (n-1)^2 / n )")
    print("-" * 20)
    print(f"Calculating for n = {n}:")
    print(f"The trace is calculated as: {c1} * ({n}+1) * ({n}-1)^2 / {n}")
    print(f"Numerator of trace = {c1} * {term1} * {term2} = {numerator_trace}")
    print(f"Denominator of trace = {denominator_trace}")
    print(f"Value of trace = {numerator_trace} / {denominator_trace} = {trace_val}")
    print(f"phi({n}) = exp({trace_val})")
    print(f"Final value: {phi_n}")

# --- User execution area ---
# Set the desired value for n (must be >= 5)
n_value = 5
calculate_phi(n_value)