import math

def calculate_periods():
    """
    Calculates the periods for the four cases based on the number of points.
    The period P is hypothesized to be sqrt(N), where N is the number of points.
    """
    # Number of points in each sequence
    n_points = [3, 4, 6, 13]
    
    # Calculate periods
    periods = [math.sqrt(n) for n in n_points]
    
    # Format the output string
    # We will format them as "sqrt(N)" for exactness, or as integer if N is a perfect square.
    period_str_parts = []
    for n in n_points:
        sqrt_n = math.sqrt(n)
        if int(sqrt_n) ** 2 == n:
            period_str_parts.append(f"{int(sqrt_n)}")
        else:
            period_str_parts.append(f"sqrt({n})")
            
    print(f"The four periods are: {', '.join(period_str_parts)}")
    # For numerical result:
    # print(f"Numerical values: {periods[0]:.4f}, {periods[1]:.4f}, {periods[2]:.4f}, {periods[3]:.4f}")

calculate_periods()
<<<sqrt(3), 2, sqrt(6), sqrt(13)>>>