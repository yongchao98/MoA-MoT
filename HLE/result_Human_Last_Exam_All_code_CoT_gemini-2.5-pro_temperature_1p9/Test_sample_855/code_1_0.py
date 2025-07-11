import math

def get_true_force():
    """Calculates the force using standard floating-point arithmetic for comparison."""
    # Planet parameters
    G = 6.6743e-11
    pi = math.pi
    r_core = 1.0e5  # m
    rho_core = 1200  # kg/m^3
    r_planet = 1.995e6 # Using average radius
    rho_shell = 300  # kg/m^3
    
    # Mass Calculation
    V_core = (4/3) * pi * r_core**3
    M_core = V_core * rho_core
    
    V_planet = (4/3) * pi * r_planet**3
    V_shell = V_planet - V_core
    M_shell = V_shell * rho_shell
    
    M_planet = M_core + M_shell
    
    # Gravity 'g'
    g = G * M_planet / r_planet**2
    
    # Force Calculation
    m_p = 50.0 # kg
    a_decel = 9.0 # m/s^2
    
    F_true = m_p * (a_decel + g)
    return F_true, g

def main():
    """
    Simulates the calculation of the landing force using Titan 5-bit architecture rules.
    """
    print("Goal: Calculate the rocket force F = m * (a + g) using Titan architecture.\n")

    # Step 1: Calculate required deceleration 'a'
    print("Step 1: Calculate the required deceleration 'a'.")
    print("The probe must decelerate from 300 m/s to 0 m/s over 5000 m.")
    print("a = v_0^2 / (2*d)")
    print("Using Titan numbers:")
    print("  v_0 = 300 -> (3/1) * 10^2")
    print("  d = 5000 -> (5/1) * 10^3")
    print("  a = ((3/1)*10^2)^2 / (2 * (5/1)*10^3) = ((9/1)*10^4) / ((10/1)*10^3)")
    print("  a = ((9/1)*10^4) / ((1/1)*10^4) = 9/1 = 9 m/s^2")
    a_val = 9

    # Step 2: Calculate the deceleration force F_decel = m * a
    print("\nStep 2: Calculate the deceleration force F_decel = m * a.")
    print("  m = 50 kg, represented as (5/1) * 10^1.")
    print("  a = 9, represented as (9/1).")
    print("  F_decel = ((5/1)*10^1) * (9/1) = (45/1) * 10^1.")
    print("The numerator 45 is larger than 31 and must be simplified.")
    print("To simplify, we renormalize the value 45 * 10^1 = 450.")
    print("  450 = 4.5 * 10^2. The fraction for 4.5 is 9/2.")
    print("  F_decel is represented as (9/2) * 10^2, which equals 450 N.")
    f_decel_val = 450
    
    # Step 3: Calculate the gravitational force F_g = m * g
    true_force, true_g = get_true_force()
    print(f"\nStep 3: Calculate the gravitational force F_g = m * g.")
    print(f"The actual gravitational acceleration 'g' on Pandora is ~{true_g:.4f} m/s^2.")
    print("To represent 'g' with Titan, we find a simple fractional approximation: 1/6.")
    print("  g ≈ 1/6 (error is less than 1%)")
    print("Now, we calculate F_g = 50 * (1/6) = 50/6.")
    print("The fraction 50/6 has a numerator > 31, so we must simplify it to 25/3.")
    print("  F_g is represented as (25/3), which equals 8.333... N.")

    # Step 4: Calculate total force by summing the components
    print("\nStep 4: Calculate the total force F_rocket = F_decel + F_g.")
    print("  F_rocket = (9/2)*10^2 + (25/3)")
    print("To add these, they need a common exponent. Converting (9/2)*10^2 to exp=0 results in a numerator of 450.")
    print("Since 450 > 31, the operation is invalid.")
    print("According to Titan rules, we must sacrifice precision by eliminating the negligible term.")
    print("Since F_g (8.33 N) << F_decel (450 N), F_g is eliminated.")
    
    f_titan = f_decel_val
    print(f"\nFinal calculated force F_titan = F_decel = {f_titan} N.")
    
    # Step 5: Calculate error
    print(f"\nStep 5: Determine the absolute error.")
    print(f"  True Force = {true_force:.3f} N")
    print(f"  Titan Calculated Force = {f_titan} N")
    error = abs(true_force - f_titan)
    print(f"  Absolute Error = |{true_force:.3f} - {f_titan}| = {error:.3f}")
    
    # Final Answer Block
    print(f"\nFinal Answer: <<<Y[{error:.3f}]>>>")

if __name__ == '__main__':
    main()