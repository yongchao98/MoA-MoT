import math

def calculate_area_contribution(n, m):
    """
    Calculates the area of the region R within a specific unit square [n, n+1) x [m, m+1).
    The area is the intersection of this square with a disk.
    """
    # Check if sqrt(n^2 + m^2) is an integer, which is a necessary condition.
    k_squared = n*n + m*m
    if k_squared < 0: return 0.0
    k_int = int(math.sqrt(k_squared))
    if k_int*k_int != k_squared:
        return 0.0
    
    k = float(k_int)
    R = k + 1.0

    def primitive_integral(x, R, C):
        """Calculates the primitive of sqrt(R^2 - t^2) - C."""
        # Clamp x to the domain of the functions
        if abs(x) > R: x = math.copysign(R, x)
        
        # Clamp asin argument
        asin_arg = x / R
        if asin_arg > 1.0: asin_arg = 1.0
        if asin_arg < -1.0: asin_arg = -1.0
        
        # Clamp sqrt argument
        sqrt_arg = R*R - x*x
        if sqrt_arg < 0: sqrt_arg = 0.0
        
        term1 = x * math.sqrt(sqrt_arg) / 2.0
        term2 = R*R * math.asin(asin_arg) / 2.0
        return term1 + term2 - C * x

    def definite_integral(a, b, R, C):
        """Calculates the definite integral of sqrt(R^2 - t^2) - C from a to b."""
        if b <= a:
            return 0.0
        return primitive_integral(b, R, C) - primitive_integral(a, R, C)

    # The area is given by integrating max(0, min(m+1, sqrt(R^2-x^2)) - m) from n to n+1.
    # We find the x-coordinates where the circle crosses the horizontal lines y=m and y=m+1.
    
    sqrt_arg_x1 = R*R - (m+1)**2
    x1 = math.sqrt(sqrt_arg_x1) if sqrt_arg_x1 > 0 else 0.0
        
    sqrt_arg_x2 = R*R - m*m
    x2 = math.sqrt(sqrt_arg_x2) if sqrt_arg_x2 > 0 else 0.0
        
    low, high = float(n), float(n+1)
    
    # The integral is split into two parts:
    # 1. A rectangular part where the integrand is 1.
    i1_low = max(low, 0.0)
    i1_high = min(high, x1)
    area1 = max(0.0, i1_high - i1_low)

    # 2. A part under the circular arc, where the integrand is sqrt(R^2-x^2) - m.
    i2_low = max(low, x1)
    i2_high = min(high, x2)
    area2 = definite_integral(i2_low, i2_high, R, float(m))
    
    return area1 + area2

if __name__ == '__main__':
    total_area = 0.0
    contributions = {}
    
    for n_val in range(6):
        for m_val in range(6):
            contrib = calculate_area_contribution(n_val, m_val)
            if contrib > 0:
                total_area += contrib
                contributions[(n_val, m_val)] = contrib
    
    A00 = contributions.get((0,0), 0)
    A10 = contributions.get((1,0), 0)
    A20 = contributions.get((2,0), 0)
    A30 = contributions.get((3,0), 0)
    A40 = contributions.get((4,0), 0)
    A50 = contributions.get((5,0), 0)
    A34 = contributions.get((3,4), 0)
    A43 = contributions.get((4,3), 0)

    print("The total area is the sum of contributions from each valid unit square.")
    print("The non-zero area contributions are:")
    print(f"A(0,0) = {A00:.4f}")
    print(f"A(1,0) = A(0,1) = {A10:.4f}")
    print(f"A(2,0) = A(0,2) = {A20:.4f}")
    print(f"A(3,0) = A(0,3) = {A30:.4f}")
    print(f"A(4,0) = A(0,4) = {A40:.4f}")
    print(f"A(5,0) = A(0,5) = {A50:.4f}")
    print(f"A(3,4) = {A34:.4f}")
    print(f"A(4,3) = {A43:.4f}")
    
    print("\nThe total area is calculated as:")
    # A00 + 2*(A10+A20+A30+A40+A50) + A34 + A43
    print(f"Area = {A00:.4f} + 2*({A10:.4f} + {A20:.4f} + {A30:.4f} + {A40:.4f} + {A50:.4f}) + {A34:.4f} + {A43:.4f}")
    
    sum_of_A_n0 = A10 + A20 + A30 + A40 + A50
    final_area = A00 + 2 * sum_of_A_n0 + A34 + A43
    
    print(f"Area = {A00:.4f} + 2*({sum_of_A_n0:.4f}) + {A34:.4f} + {A43:.4f}")
    print(f"Area = {A00:.4f} + {2*sum_of_A_n0:.4f} + {A34:.4f} + {A43:.4f}")
    print(f"Area = {final_area:.2f}")
    print(f"<<<{final_area:.2f}>>>")