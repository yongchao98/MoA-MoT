import math

def solve_parrot_math():
    """
    This function instructs a parrot to estimate the mass of a rock.
    It follows these steps:
    1. Define the given values as fractions with integers <= 10.
    2. Choose an approximation for Pi that satisfies the error and integer constraints.
    3. Construct and print the final calculation equation.
    """
    
    # Step 1: Define values as fractions.
    # Density (rho) = 0.9 = 9/10 kg/cm^3
    density_num, density_den = 9, 10
    
    # Radius (r) = 0.5 = 1/2 cm
    radius_num, radius_den = 1, 2
    
    # Constant from volume formula (4/3)
    vol_const_num, vol_const_den = 4, 3
    
    # Step 2: Choose approximation for Pi.
    # We choose pi ~ 3/1 because it uses smaller integers than 10/3
    # and keeps the error within 10%.
    # Error for pi ~ 3/1 is approx. 4.5%
    pi_approx_num, pi_approx_den = 3, 1
    
    # Radius cubed
    r_cubed_num = radius_num ** 3
    r_cubed_den = radius_den ** 3

    # All integers used in the calculation:
    # From density: 9, 10
    # From volume constant: 4, 3
    # From pi approximation: 3, 1
    # From radius: 1, 2
    # The set of unique integers is {1, 2, 3, 4, 9, 10}. All are <= 10.
    
    # Step 3: Print the final calculation for the parrot.
    print("Yes, the parrot can estimate the mass.")
    print("Here is the calculation using fractions with integers no larger than 10:")
    print("Mass = Density * (4/3) * pi * (Radius)^3")
    print("\nSubstitute the fractional values:")
    
    # Print the equation with each number, as requested.
    print(f"Mass ≈ {density_num}/{density_den} * {vol_const_num}/{vol_const_den} * {pi_approx_num}/{pi_approx_den} * ({radius_num}/{radius_den})^3")
    
    # Calculate the result to show the parrot the final answer
    final_num = density_num * vol_const_num * pi_approx_num * r_cubed_num
    final_den = density_den * vol_const_den * pi_approx_den * r_cubed_den
    
    # Simplify the fraction by finding the greatest common divisor
    common_divisor = math.gcd(final_num, final_den)
    simple_num = final_num // common_divisor
    simple_den = final_den // common_divisor
    
    print(f"\nThis simplifies to the final estimated mass of {simple_num}/{simple_den} kg.")

solve_parrot_math()