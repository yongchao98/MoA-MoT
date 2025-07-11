import math

def is_inside(a, n):
    """
    Checks if the pole z = a + 2*pi*i*n is inside the curve Gamma.
    The point (x,y) is (a, 2*pi*n).
    """
    y = 2 * math.pi * n
    u = a + y
    v = a - y
    
    # The equation for the curve is G(u,v) = 0
    # The inside is where G(u,v) < 0
    val = 3 * (u**2 - 400)**3 - 20 * u**3 * v**2 + 1200 * v**4
    return val < 0

def solve():
    """
    Calculates the contour integral by counting the poles inside the curve.
    """
    pole_count = 0
    # The curve is bounded, let's check a reasonable range for a and n.
    # From manual analysis, intercepts are around +/- 13.5, so |y| < 14.
    # |2*pi*n| < 14  => |n| < 14/(2*pi) approx 2.23. So n=-2, -1, 0, 1, 2.
    # |x| = |a| < 14. Let's check a wider range for safety.
    
    a_min = -2024
    a_max = 2024
    n_max = 2 # Based on |y| < 14

    for n in range(-n_max, n_max + 1):
        for a in range(a_min, a_max + 1):
            if is_inside(a, n):
                pole_count += 1
                
    # We can optimize the search for a, as they form contiguous intervals for each n.
    # For a given n, find the start and end of the integer range for 'a'.
    
    pole_count_optimized = 0
    for n in range(-n_max, n_max + 1):
        # Find the first 'a' from the left that is inside
        start_a = None
        for a in range(-30, 31): # Search in a reasonable range around 0
            if is_inside(a, n):
                start_a = a
                break
        
        # Find the first 'a' from the right that is inside
        end_a = None
        for a in range(30, -31, -1):
            if is_inside(a, n):
                end_a = a
                break
        
        if start_a is not None and end_a is not None:
            # Add the number of integers in the range [start_a, end_a]
            # that are also within the function's definition range [-2024, 2024]
            actual_start = max(start_a, a_min)
            actual_end = min(end_a, a_max)
            if actual_end >= actual_start:
                 pole_count_optimized += (actual_end - actual_start + 1)
                 
    # The integral is 2 * pi * i * pole_count
    # The question asks for the value of the integral.
    # We print the components of the final equation.
    
    # Let's run the optimized count to be sure.
    n_counts = {}
    for n in range(-n_max, n_max + 1):
        count = 0
        for a in range(a_min, a_max + 1):
             if is_inside(a,n):
                 count += 1
        n_counts[n] = count
    
    final_pole_count = sum(n_counts.values())

    print(f"The total number of poles inside the contour Gamma is {final_pole_count}.")
    print(f"The value of the integral is 2 * pi * i * {final_pole_count}.")
    #The result is an imaginary number.
    print(f"The final equation for the integral is 2 * 3.14159... * i * {final_pole_count} = {2*final_pole_count} * pi * i")
    print(f"The coefficient of pi*i is {2*final_pole_count}.")

solve()