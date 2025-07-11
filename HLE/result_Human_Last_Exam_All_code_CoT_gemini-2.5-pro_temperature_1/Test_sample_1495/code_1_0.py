import random

def solve_triangle_probability():
    """
    This function calculates the probability that a random point P in triangle ABC
    is inside triangle XYZ, where XYZ is formed by three random cevians.

    The probability is the expected value of the area ratio Area(XYZ)/Area(ABC).
    We use the Monte Carlo method to estimate this expected value.
    """
    
    num_samples = 10_000_000  # Number of random trials for the simulation
    total_area_ratio_sum = 0.0

    for _ in range(num_samples):
        # Step 1: Uniformly randomly choose points D, E, F on the sides
        # by generating u, v, w from a uniform distribution in [0, 1].
        u = random.random()
        v = random.random()
        w = random.random()

        # Step 2: Calculate the area ratio using the formula derived from Routh's Theorem.
        # The formula is: (uvw - (1-u)(1-v)(1-w))^2 / ((1-u+uv)(1-v+vw)(1-w+wu))
        
        # Numerator of the equation
        term1 = u * v * w
        term2 = (1 - u) * (1 - v) * (1 - w)
        numerator = (term1 - term2)**2

        # Denominator of the equation
        den_factor1 = 1 - u + u * v
        den_factor2 = 1 - v + v * w
        den_factor3 = 1 - w + w * u
        
        # This check is for numerical stability, though division by zero is highly unlikely
        # with random floating-point numbers.
        if den_factor1 == 0 or den_factor2 == 0 or den_factor3 == 0:
            continue
            
        denominator = den_factor1 * den_factor2 * den_factor3
        
        area_ratio = numerator / denominator
        total_area_ratio_sum += area_ratio

    # Step 3: The probability is the average of all calculated ratios.
    probability = total_area_ratio_sum / num_samples

    print(f"The probability that point P is in triangle XYZ is estimated to be: {probability}")

solve_triangle_probability()