import random

def get_area_ratio(u, v, w):
    """
    Calculates the ratio of Area(XYZ) / Area(ABC) based on Routh's Theorem.
    
    The inputs u, v, w are the fractional positions of points D, E, F on
    their respective sides of triangle ABC.
    u = BD / BC
    v = CE / CA
    w = AF / AB
    
    The formula for the ratio is:
    ((1-u)(1-v)(1-w) - u*v*w)**2 / ((1-u+u*v)*(1-v+v*w)*(1-w+w*u))
    """
    
    # Numerator calculation
    term1 = (1.0 - u) * (1.0 - v) * (1.0 - w)
    term2 = u * v * w
    numerator = (term1 - term2)**2
    
    # Denominator calculation
    den1 = 1.0 - u + u * v
    den2 = 1.0 - v + v * w
    den3 = 1.0 - w + w * u
    denominator = den1 * den2 * den3
    
    # Handle the rare case where the denominator is zero. This corresponds to a
    # degenerate inner triangle with zero area.
    if denominator == 0.0:
        return 0.0
        
    return numerator / denominator

def main():
    """
    Performs a Monte Carlo simulation to estimate the probability.
    """
    num_samples = 2000000  # Number of random configurations to simulate
    total_ratio_sum = 0.0
    
    for _ in range(num_samples):
        # Generate three independent random numbers from a uniform distribution [0, 1]
        u = random.random()
        v = random.random()
        w = random.random()
        
        # Calculate the area ratio for this configuration and add it to the total
        total_ratio_sum += get_area_ratio(u, v, w)
        
    # The estimated probability is the average of all calculated ratios
    estimated_probability = total_ratio_sum / num_samples
    
    # Based on the simulation, the final equation is that the expected
    # value of the area ratio is approximately the calculated number.
    # The problem asks to output the numbers in the final equation. 
    # Here, the final equation is: Probability = [calculated_value]
    
    print("The probability that point P is in triangle XYZ is estimated to be:")
    print(estimated_probability)

if __name__ == "__main__":
    main()