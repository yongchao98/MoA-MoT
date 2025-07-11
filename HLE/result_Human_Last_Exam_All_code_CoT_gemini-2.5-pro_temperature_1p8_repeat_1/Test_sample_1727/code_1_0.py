import math

def calculate_chi(kh):
    """
    Calculates the stable amplitude factor chi based on the dimensionless parameter kh.

    Args:
        kh (float): The dimensionless parameter, product of wavenumber k and thickness h.

    Returns:
        float: The value of the amplitude factor chi.
        Or None if the denominator is zero.
    """
    
    # Calculate sinh(kh) and cosh(kh)
    try:
        sinh_kh = math.sinh(kh)
        cosh_kh = math.cosh(kh)
    except OverflowError:
        print(f"Error: kh value {kh} is too large and caused a math overflow.")
        return None

    # Numerator of the expression for chi
    numerator = kh
    
    # Denominator of the expression for chi
    denominator = kh * sinh_kh - cosh_kh
    
    # Avoid division by zero, which happens at resonance (kh * tanh(kh) = 1)
    if abs(denominator) < 1e-9:
        print(f"The value kh = {kh} is at or near a resonance point. Chi is infinite.")
        return None
        
    # Calculate chi
    chi = numerator / denominator

    # Print the final equation with all numbers
    print("The formula for the amplitude factor is: chi = kh / (kh * sinh(kh) - cosh(kh))")
    print(f"\nFor the case where kh = {kh}:")
    
    # Print each part of the calculation
    print(f"Numerator = kh = {numerator}")
    
    print(f"Denominator = ({kh} * sinh({kh})) - cosh({kh})")
    print(f"            = ({kh} * {sinh_kh:.6f}) - {cosh_kh:.6f}")
    print(f"            = {kh * sinh_kh:.6f} - {cosh_kh:.6f}")
    print(f"            = {denominator:.6f}")

    print(f"\nResult: chi = {numerator} / {denominator:.6f} = {chi:.6f}")
    
    return chi

# --- Main execution ---
# We will evaluate chi for a representative value, e.g., kh = 2.0
# This value is chosen for demonstration purposes.
kh_value = 2.0
final_chi = calculate_chi(kh_value)

# The problem is solved analytically, resulting in a formula for chi.
# The numerical value depends on the specific value of kh.
# For kh=2.0, the value of chi is calculated above.
# The final result is the expression for chi as a function of kh.
# Here we provide the numerical result for the chosen representative value.
if final_chi is not None:
    # The problem asks for the stable amplitude factor. Here we provide its numerical value
    # for our chosen case, kh=2.0. The formula derived is general.
    # We will output the numerical value in the requested format.
    print(f"\n<<<Amplitude factor for kh={kh_value} is {final_chi:.4f}>>>")
