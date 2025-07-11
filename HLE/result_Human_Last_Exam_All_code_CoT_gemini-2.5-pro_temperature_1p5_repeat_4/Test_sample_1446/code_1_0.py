def calculate_mean_field_nu():
    """
    Calculates the critical exponent ν within the mean-field approximation
    of the G₄ (φ⁴) theoretical framework.

    In the Ginzburg-Landau-Wilson (φ⁴) theory of critical phenomena, the critical
    exponents take on their classical, mean-field values for spatial dimensions 'd'
    at or above the upper critical dimension (d_c = 4).

    The critical exponent ν governs the scaling of the correlation length (ξ)
    near the critical temperature (T_c): ξ ~ |T - T_c|^(-ν).

    The mean-field value for ν is derived from the Landau free energy expansion
    and is precisely 1/2.
    """
    
    # Define the numbers in the equation for ν
    numerator = 1
    denominator = 2
    
    # Calculate the value of the critical exponent ν
    nu = numerator / denominator
    
    # Print the context and the equation step-by-step
    print("Within the mean-field approximation of the G₄ (φ⁴) theory (valid for dimensions d >= 4), the critical exponent ν is determined by the following equation:")
    print(f"ν = {numerator} / {denominator}")
    
    # Print the final precise value
    print(f"The precise value is: {nu}")

if __name__ == "__main__":
    calculate_mean_field_nu()