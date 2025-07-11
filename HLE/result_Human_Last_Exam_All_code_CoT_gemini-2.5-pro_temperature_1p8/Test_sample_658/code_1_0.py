import math

def calculate_asymptotic_degree(B, delta):
    """
    Calculates the asymptotic value of d_{B,delta}.

    Args:
        B (float): The upper bound of the interval [0, B], B >= 1.
        delta (float): The approximation error, 0 < delta < 1.
    
    Returns:
        float: The asymptotic value A(B, delta).
    """
    if B < 1 or delta <= 0 or delta >= 1:
        raise ValueError("Inputs are out of the valid range.")

    L = math.log(1.0 / delta)

    # To avoid issues with log(x) where x can be non-positive, 
    # and to capture the behavior correctly, we use log(L/B + C) with C > 1.
    # We choose C=2 for robustness.
    
    # Asymptotic behavior for L/B -> 0
    term1 = math.sqrt(B * L)
    
    # Asymptotic behavior for L/B -> infinity
    denominator_log_arg = L / B + 2.0
    term2 = L / math.log(denominator_log_arg)
    
    asymptotic_value = term1 + term2
    
    print(f"Given B = {B} and delta = {delta}")
    print(f"L = log(1/delta) = {L}")
    print(f"Term 1 (sqrt(B*L)): {term1}")
    print(f"Term 2 (L / log(L/B + 2)): {term2}")
    print(f"Asymptotic value d_{{B,delta}} = Theta(A(B,delta)) where A(B,delta) is: {asymptotic_value}")

if __name__ == '__main__':
    # Example values, can be changed by the user.
    # Regime 1: B >> L
    # B_val_1 = 1000
    # delta_val_1 = 0.01 
    # calculate_asymptotic_degree(B_val_1, delta_val_1)
    
    # Regime 2: L >> B
    # B_val_2 = 2
    # delta_val_2 = 1e-100
    # calculate_asymptotic_degree(B_val_2, delta_val_2)

    # Regime 3: L and B are comparable
    B_val_3 = 10
    delta_val_3 = 1e-4
    calculate_asymptotic_degree(B_val_3, delta_val_3)
