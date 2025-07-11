def calculate_ratio():
    """
    Calculates the ratio R based on the one-loop counter-term coefficients.
    
    The counter-terms are expressed in units of a common factor C = g^2 / (32*pi^2*epsilon),
    which cancels out in the final ratio.
    """
    # Coefficients of the common factor C for each counter-term
    delta_Zx_coeff = 1
    delta_Zmx_coeff = -3
    delta_Zg_coeff = 2

    # The ratio R is the ratio of these coefficients
    numerator = delta_Zx_coeff
    denominator = delta_Zg_coeff + delta_Zmx_coeff
    
    R = numerator / denominator
    
    print("The counter-term coefficients are proportional to a common factor C:")
    print(f"delta_Zx = {delta_Zx_coeff} * C")
    print(f"delta_Zg = {delta_Zg_coeff} * C")
    print(f"delta_Zmx = {delta_Zmx_coeff} * C")
    print("\nThe ratio R is calculated as:")
    print(f"R = delta_Zx / (delta_Zg + delta_Zmx)")
    print(f"R = ({delta_Zx_coeff}*C) / (({delta_Zg_coeff}*C) + ({delta_Zmx_coeff}*C))")
    print(f"R = {delta_Zx_coeff} / ({delta_Zg_coeff} + {delta_Zmx_coeff})")
    print(f"R = {numerator} / ({denominator})")
    print(f"R = {R}")

if __name__ == "__main__":
    calculate_ratio()