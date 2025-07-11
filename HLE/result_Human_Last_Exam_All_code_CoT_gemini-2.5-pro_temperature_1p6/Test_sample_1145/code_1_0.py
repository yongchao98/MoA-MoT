import math

def calculate_probability():
    """
    Calculates the probability of finding a particle in a 1D box
    in a given interval for n=2.
    """
    # Given parameters
    n = 2
    x1_norm = 0.495  # x1/a
    x2_norm = 0.505  # x2/a

    # The probability is the definite integral of the probability density:
    # P = Integral from x1 to x2 of (2/a) * sin^2(n*pi*x/a) dx
    # The antiderivative (evaluated in terms of x/a) is:
    # F(x/a) = (x/a) - (1 / (2*n*pi)) * sin(2*n*pi*x/a)
    # The definite integral is F(x2/a) - F(x1/a).

    # Let's calculate the terms of the final equation:
    # P = (x2/a - x1/a) - (1/(2*n*pi)) * [sin(2*n*pi*x2/a) - sin(2*n*pi*x1/a)]
    
    term1 = x2_norm - x1_norm
    
    factor = 1 / (2 * n * math.pi)
    
    arg1 = 2 * n * math.pi * x1_norm
    sin_term1 = math.sin(arg1)
    
    arg2 = 2 * n * math.pi * x2_norm
    sin_term2 = math.sin(arg2)
    
    sin_diff = sin_term2 - sin_term1
    
    term2 = factor * sin_diff

    probability = term1 - term2

    # --- Output the results ---
    print("The probability P for n=2 in the interval [0.495a, 0.505a] is given by the equation:")
    print("P = (x2/a - x1/a) - (1/(2*n*pi)) * (sin(2*n*pi*x2/a) - sin(2*n*pi*x1/a))")
    print("\nSubstituting the values:")
    
    # We show the numeric components of the equation
    # P = (0.505 - 0.495) - (1/(4*pi)) * (sin(4*pi*0.505) - sin(4*pi*0.495))
    # Note: 4*pi*0.505 = 2.02*pi and 4*pi*0.495 = 1.98*pi
    print(f"P = ({x2_norm} - {x1_norm}) - (1/({2*n}*pi)) * (sin({2*n}*pi*{x2_norm}) - sin({2*n}*pi*{x1_norm}))")
    print(f"P = {term1} - {factor:.6f} * ({sin_term2:.6f} - ({sin_term1:.6f}))")
    print(f"P = {term1} - {factor:.6f} * ({sin_diff:.6f})")
    print(f"P = {term1} - {term2:.9f}")
    
    print(f"\nFinal Probability: {probability}")
    
if __name__ == '__main__':
    calculate_probability()
<<<6.579544552216591e-06>>>