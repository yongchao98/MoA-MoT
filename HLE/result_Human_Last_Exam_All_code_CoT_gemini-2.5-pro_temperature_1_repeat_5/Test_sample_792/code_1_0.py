import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def multiply_fractions(f1_num, f1_den, f2_num, f2_den):
    """
    Multiplies two fractions according to Titan's rules.
    Raises a ValueError if the calculation is invalid.
    """
    # Pre-simplify by dividing out common factors
    g1 = gcd(f1_num, f2_den)
    g2 = gcd(f2_num, f1_den)
    
    num1_s = f1_num // g1
    den2_s = f2_den // g1
    
    num2_s = f2_num // g2
    den1_s = f1_den // g2

    # Check if intermediate multiplications exceed 31
    # This check is implicitly handled by the final check below in this case
    
    # Calculate final numerator and denominator
    res_num = num1_s * num2_s
    res_den = den1_s * den2_s
    
    # Check if the final result respects the 5-bit constraint
    if not (0 <= res_num <= 31 and 0 <= res_den <= 31):
        raise ValueError(f"Result {res_num}/{res_den} is out of 5-bit range.")
        
    return res_num, res_den

def solve():
    """
    Solves the Curious Monkey and Lion problem using Titan's computational rules.
    """
    # 1. Physics Derivation (performed offline)
    # Force F = (d * m * g) / h
    # Given d=20m, h=10m, g≈10m/s^2, r=0.5cm, ρ=0.9kg/cm^3
    # Mass m = Volume * density = (4/3 * pi * r^3) * ρ
    # m = (4/3 * pi * (0.5cm)^3) * (0.9 kg/cm^3)
    # m = (4/3 * pi * 0.125) * 0.9 kg = (pi/6) * 0.9 kg = 0.15 * pi kg
    # F = (20m * (0.15*pi kg) * 10m/s^2) / 10m
    # F = 20 * 0.15 * pi N = 3 * pi N
    # The task for Titan is to calculate F = 3 * pi.

    # 2. Find the best approximation for pi
    # We need to find n/d ≈ pi such that 3/1 * n/d is a valid Titan operation.
    # The rule for a/b * c/d requires (a/gcd(a,d)) * (c/gcd(c,b)) <= 31.
    # For 3/1 * n/d, this means (3/gcd(3,d)) * n <= 31.
    # To allow for a larger n (and thus better accuracy), gcd(3,d) should be 3.
    # This means d must be a multiple of 3.
    # Constraint becomes (3/3) * n <= 31, so n <= 31.
    # We search for the best fraction n/d ≈ pi with d being a multiple of 3.
    # - pi ≈ 19/6 ≈ 3.1667. Calculation: 3/1 * 19/6 -> 19/2 = 9.5
    # - pi ≈ 22/7 (Invalid, 3*22=66>31)
    # - pi ≈ 28/9 ≈ 3.1111. Calculation: 3/1 * 28/9 -> 28/3 ≈ 9.333
    # The approximation 19/6 gives a result of 9.5, with error |9.5 - 3*pi| ≈ 0.075
    # The approximation 28/9 gives a result of 9.333, with error |9.333 - 3*pi| ≈ 0.091
    # The best approximation is pi ≈ 19/6.
    
    pi_approx_num, pi_approx_den = 19, 6
    three_num, three_den = 3, 1

    # 3. Perform the Titan calculation
    try:
        force_num, force_den = multiply_fractions(
            three_num, three_den, 
            pi_approx_num, pi_approx_den
        )

        # 4. Report the results
        print("Yes, the force can be calculated.")
        print("The calculation performed is F = 3 * pi.")
        print(f"Using the fractional approximation pi ≈ {pi_approx_num}/{pi_approx_den}:")
        
        # Print the final equation with all its numbers
        print("\nFinal Titan Equation:")
        print(f"{force_num} / {force_den} = {three_num} / {three_den} * {pi_approx_num} / {pi_approx_den}")
        
        calculated_force = force_num / force_den
        true_force = 3 * math.pi
        absolute_error = abs(calculated_force - true_force)
        
        print(f"\nCalculated Force: {calculated_force:.4f} N")
        print(f"True Force: {true_force:.4f} N")
        print(f"Smallest Absolute Error (e): {absolute_error:.6f}")
        
        # This is the final answer for the submission system
        answer = f"Y[{absolute_error:.6f}]"
        # print(f"\nAnswer string: <<< {answer} >>>")
        print(f"<<<{answer}>>>")


    except ValueError as e:
        print(f"The calculation is not possible. Reason: {e}")
        print("<<<N0>>>")

solve()
