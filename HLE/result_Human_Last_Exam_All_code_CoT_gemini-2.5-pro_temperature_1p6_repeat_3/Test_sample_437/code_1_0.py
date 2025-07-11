import math

def evaluate_1s_integral():
    """
    Evaluates the integral <phi_1s | 1/r | phi_1s> for a 1s Slater orbital.
    
    The script explains and calculates the value step-by-step.
    """
    
    # Let's use a specific value for the orbital exponent zeta for demonstration.
    # zeta = 1.0 is the exact value for the Hydrogen atom's 1s orbital.
    zeta = 1.0
    
    print("--- Evaluating the integral <phi_1s | 1/r | phi_1s> ---")
    print("\nThis integral represents the expectation value of the potential energy operator 1/r.")
    print(f"We will use a 1s Slater-type orbital with exponent zeta = {zeta:.2f}.\n")
    
    # --- The Integral Setup ---
    print("The integral in spherical coordinates is:")
    print("I = Integral( [phi_1s(r)]^2 * (1/r) * r^2 * sin(theta) dr dtheta dphi )")
    print("where phi_1s(r) = N * exp(-zeta*r) and N^2 = zeta^3 / pi.\n")
    
    print("This can be broken down into three parts:")
    print("I = (Constant Part) * (Angular Integral) * (Radial Integral)\n")
    
    # --- Step 1: Constant Part ---
    constant_part = (zeta**3) / math.pi
    print(f"Step 1: The Constant Part (from N^2)")
    print(f"N^2 = zeta^3 / pi = {zeta**3:.2f} / {math.pi:.4f} = {constant_part:.4f}")
    
    # --- Step 2: Angular Integral ---
    angular_integral = 4 * math.pi
    print("\nStep 2: The Angular Integral")
    print("For any spherically symmetric function, the integral over the solid angle is 4*pi.")
    print(f"Angular Integral = 4 * pi = {angular_integral:.4f}")
    
    # --- Step 3: Radial Integral ---
    # The radial integrand is exp(-2*zeta*r) * (1/r) * r^2 = r * exp(-2*zeta*r)
    # The integral is of the standard form Integral(x^n * exp(-ax) dx) from 0 to infinity,
    # which evaluates to n! / a^(n+1).
    # Here, n=1 and a = 2*zeta.
    n = 1
    a = 2 * zeta
    radial_integral = math.factorial(n) / (a**(n + 1))
    print("\nStep 3: The Radial Integral")
    print("The radial part is Integral( r * exp(-2*zeta*r) dr ) from 0 to infinity.")
    print(f"Using the standard formula with n={n} and a=2*zeta={a:.2f}, the result is:")
    print(f"Radial Integral = n! / a^(n+1) = {math.factorial(n)} / ({a:.2f})^({n+1}) = {radial_integral:.4f}")
    
    # --- Step 4: Final Calculation ---
    final_result = constant_part * angular_integral * radial_integral
    print("\n--- Final Calculation ---")
    print("Now, we multiply the parts together.")
    print("Result = (Constant Part) * (Angular Integral) * (Radial Integral)")
    
    # Print the final equation with the numbers
    print("\nFinal equation with numerical values:")
    print(f"Result = {constant_part:.4f} * {angular_integral:.4f} * {radial_integral:.4f}")
    print(f"Result = {final_result:.4f}\n")
    
    print("Analytically, the result is (zeta^3 / pi) * (4 * pi) * (1 / (4 * zeta^2)) = zeta.")
    print(f"For our chosen zeta = {zeta:.2f}, the analytical result is {zeta:.4f}, which matches our calculation.")

if __name__ == '__main__':
    evaluate_1s_integral()
    # The final result is zeta, which for our example is 1.0.
    # The problem asks for the evaluation of the expression. For zeta=1.0, the value is 1.0.
    print("\n<<<1.0>>>")
