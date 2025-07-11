def calculate_ratio():
    """
    Calculates the ratio R of one-loop counter-terms in Yukawa theory.

    The counter-terms at one-loop in the MS-bar scheme are:
    delta_Zx = (g^2 / (32 * pi^2 * epsilon))
    delta_Zmx = -3 * (g^2 / (32 * pi^2 * epsilon))
    delta_Zg = (g^2 / (32 * pi^2 * epsilon))
    
    We can define a common factor C = g^2 / (32 * pi^2 * epsilon) and express
    the counter-terms as integer multiples of C. The factor C will cancel out
    in the ratio R.
    """
    
    # Coefficients of the common factor C for each counter-term
    coeff_dZx = 1
    coeff_dZmx = -3
    coeff_dZg = 1
    
    # The numerator of R is delta_Zx
    numerator = coeff_dZx
    
    # The denominator of R is delta_Zg + delta_Zmx
    denominator = coeff_dZg + coeff_dZmx
    
    # Calculate the ratio R
    if denominator == 0:
        R = "undefined (division by zero)"
    else:
        R = numerator / denominator

    # Print the equation with the numerical coefficients
    print(f"The ratio R is calculated as:")
    print(f"R = delta_Zx / (delta_Zg + delta_Zmx)")
    print(f"Using the coefficients of the common factor C = g^2/(32*pi^2*epsilon):")
    print(f"R = {numerator} / ({coeff_dZg} + ({coeff_dZmx}))")
    
    # Print the final calculation and result
    print(f"R = {numerator} / {denominator} = {R}")

calculate_ratio()