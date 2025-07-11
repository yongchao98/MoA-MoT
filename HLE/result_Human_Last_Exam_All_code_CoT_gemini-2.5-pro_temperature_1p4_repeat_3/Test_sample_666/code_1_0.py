import math

def count_poles_inside_curve():
    """
    Counts the number of poles of the function f(z) that lie inside the curve Gamma.

    The curve is defined by the equation P(x,y) = 0, which can be expressed in
    u=x+y, v=x-y coordinates as Q(u,v) = 3(u^2-400)^3 - 20u^3v^2 + 1200v^4 = 0.
    A point is inside the curve if Q(u,v) < 0.

    The poles are located at z = a + 2*pi*i*k, where 'a' and 'k' are integers.
    This corresponds to x=a and y=2*pi*k.
    In the transformed coordinates, u = a + 2*pi*k and v = a - 2*pi*k.

    Analysis shows the curve is bounded, so we only need to check a small range
    of 'a' and 'k' values.
    """
    count = 0
    pi = math.pi
    
    # Range for a is given as [-2024, 2024], but analysis shows we only need to check |a| < 22.
    # We use a slightly larger range for safety.
    a_range = range(-25, 26)
    
    # Analysis shows |k| must be <= 3. We use a slightly larger range for safety.
    k_range = range(-4, 5)

    for a in a_range:
        for k in k_range:
            x = float(a)
            y = 2.0 * pi * k
            
            u = x + y
            v = x - y
            
            # Evaluate the polynomial Q(u,v)
            try:
                u2 = u * u
                v2 = v * v
                
                term1 = 3.0 * (u2 - 400.0)**3
                term2 = -20.0 * u**3 * v2
                term3 = 1200.0 * v2**2
                
                q_value = term1 + term2 + term3
            except OverflowError:
                # For large u,v this value is large and positive.
                q_value = 1.0e30 # Represents a large positive number
            
            if q_value < 0:
                count += 1
    
    return count

def main():
    """
    Calculates the contour integral and prints the result.
    """
    # The integral is 2 * pi * i * N, where N is the number of enclosed poles.
    # The residue at each pole is 1.
    num_poles = count_poles_inside_curve()
    
    # The contour integral is 2 * pi * i * num_poles
    coefficient = 2 * num_poles
    
    # The problem asks for the contour integral. We present it symbolically.
    # The original polynomial coefficients are:
    c1 = 3
    c2 = 18
    c3 = 45
    c4 = 60
    # ... and so on. The instruction "output each number in the final equation"
    # is interpreted as showing the final calculation clearly.
    
    print(f"The number of enclosed poles is: {num_poles}")
    print(f"The value of the contour integral is 2 * {num_poles} * pi * i.")
    print(f"Final answer: {coefficient}*pi*i")

if __name__ == "__main__":
    main()
