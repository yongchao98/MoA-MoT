import math

def find_critical_n():
    """
    Calculates the critical Hill coefficient (n_c) for the onset of oscillations.
    """
    # System parameters
    alpha = 100.0
    beta = 20.0

    print("--- Calculating the critical Hill coefficient n_c ---")
    print(f"Given parameters: alpha = {alpha}, beta = {beta}")
    print("\nStep 1: Calculate the steady-state concentration x_s at the bifurcation point.")
    print("Formula: x_s = sqrt(beta) - 1")
    sqrt_beta = math.sqrt(beta)
    x_s = sqrt_beta - 1
    print(f"sqrt({beta}) = {sqrt_beta}")
    print(f"x_s = {sqrt_beta} - 1 = {x_s}")

    print("\nStep 2: Calculate n_c using the derived formula.")
    print("Formula: n_c = log(alpha / (2*sqrt(beta) - 1) - 1) / log(sqrt(beta) - 1)")
    
    # Calculate terms for the formula
    numerator_arg_term1 = alpha / (2 * sqrt_beta - 1)
    numerator_arg = numerator_arg_term1 - 1
    denominator_arg = x_s
    
    # Check if oscillations are possible
    if numerator_arg <= 0:
        print("\nOscillations are not possible for n > 0 with the given parameters.")
        return

    # Perform the calculation and print intermediate steps
    print("\n--- Calculation Breakdown ---")
    print(f"Value of the numerator argument for the log function:")
    print(f"alpha / (2*sqrt(beta) - 1) - 1 = {alpha} / (2*{sqrt_beta} - 1) - 1 = {numerator_arg_term1} - 1 = {numerator_arg}")

    print(f"\nValue of the denominator argument for the log function:")
    print(f"sqrt(beta) - 1 = {denominator_arg}")
    
    log_numerator = math.log(numerator_arg)
    log_denominator = math.log(denominator_arg)
    n_c = log_numerator / log_denominator

    print(f"\nFinal calculation:")
    print(f"n_c = log({numerator_arg}) / log({denominator_arg})")
    print(f"n_c = {log_numerator} / {log_denominator}")
    print(f"n_c = {n_c}")
    
    print("\n--- Conclusion ---")
    print(f"The system exhibits oscillations for values of the Hill coefficient n > {n_c}.")

find_critical_n()