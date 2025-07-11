import math

def get_asymptotic_d(B, delta):
    """
    Calculates the asymptotic value of the minimum degree d_B,delta
    for approximating e^-x on [0, B] with error delta.
    
    Args:
        B (float): The upper bound of the interval [0, B]. Must be >= 1.
        delta (float): The maximum allowed error. Must be in (0, 1).
        
    Returns:
        float: The asymptotic value A(B, delta).
    """
    if B < 1 or not (0 < delta < 1):
        raise ValueError("B must be >= 1 and delta must be in (0, 1).")
        
    L = math.log(1 / delta)
    
    # First term of the asymptotic formula
    term1 = math.sqrt(B * L)
    
    # Second term of the asymptotic formula
    # We use log(e + L/B) to handle the case where L/B might be small
    # and to provide a smooth transition between regimes.
    term2 = L / math.log(math.e + L / B)
    
    result = term1 + term2
    
    # Output the components of the final equation as requested
    print(f"For B = {B} and delta = {delta}:")
    print(f"L = log(1/delta) = {L:.4f}")
    print(f"Term 1 (sqrt(B*L)): {term1:.4f}")
    print(f"Term 2 (L / log(e + L/B)): {term2:.4f}")
    print(f"Asymptotic value d_{B,delta} ≈ {term1:.4f} + {term2:.4f} = {result:.4f}")
    
    return result

if __name__ == '__main__':
    # Example values for B and delta
    B_val = 100
    delta_val = 0.01
    
    # Calculate and print the asymptotic value
    asymptotic_d = get_asymptotic_d(B_val, delta_val)