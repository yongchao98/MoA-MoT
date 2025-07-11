import math

# Titan Computer Simulation for the Curious Monkey Problem

def main():
    """
    Solves the physics problem using the constraints of the Titan superconducting computer.
    """
    print("### Titan Computer Calculation Plan ###")
    print("1. Establish the physics formula for the required force F.")
    print("2. Represent all physical constants and parameters as fractions of 5-bit integers (0-31).")
    print("3. Perform the calculation using fractional arithmetic, ensuring no intermediate numerator or denominator exceeds 31.")
    print("4. Output the final equation, the resulting force as a fraction, and the calculated error.")
    print("-" * 30)

    # --- Step 1: Physics Formula ---
    print("\n### Step 1: Deriving the Physics Formula ###")
    print("The time 't' for the rock to fall 10m is given by h = (1/2)gt^2.")
    print("The horizontal distance 'd' covered is given by d = (1/2)a_x * t^2, where a_x = F/m.")
    print("Combining these, we get d = (1/2)(F/m)(2h/g) = F*h/(m*g).")
    print("Solving for the force F, we get: F = (d * m * g) / h.")
    print("The mass 'm' of the spherical rock is m = density * volume = rho * (4/3) * pi * r^3.")
    print("Substituting m, the full formula is: F = (d/h) * g * rho * (4/3) * pi * r^3.")
    
    # --- Step 2: Titan Fraction Representation ---
    print("\n### Step 2: Approximating Constants as Titan Fractions ###")
    
    # Define system parameters in SI units
    d_si = 20.0  # meters
    h_si = 10.0  # meters
    g_si = 9.8   # m/s^2
    r_si = 0.005 # meters (0.5 cm)
    rho_si = 900000.0 # kg/m^3 (0.9 kg/cm^3)

    # Values for the formula components
    dist_ratio_val = d_si / h_si
    four_thirds_val = 4.0 / 3.0
    mass_factor_val = rho_si * (r_si**3) # This is rho * r^3
    pi_val = math.pi
    
    print(f"The core calculation is F = (d/h) * (4/3) * g * pi * (rho * r^3)")

    # The values for r (0.005) and rho (900,000) are outside the representable range of Titan.
    # However, we can calculate their product first and then approximate it as a Titan fraction.
    print(f"\nThe term (rho * r^3) evaluates to {mass_factor_val:.4f} kg.")
    
    # To solve this problem, we must choose fractional approximations that allow for chained
    # multiplication without intermediate overflow. This requires careful selection to enable
    # cross-cancellation at each step. After analysis, the following set of approximations is chosen.
    # While not individually the most accurate, they allow the calculation to complete successfully.

    dist_ratio_frac = (2, 1)  # d/h = 20/10 = 2/1 (Exact)
    four_thirds_frac = (4, 3) # (Exact)
    
    # Approximations chosen for computability:
    mass_factor_frac = (1, 8) # Approx for 0.1125. (0.125, 11% error)
    pi_frac = (3, 1)          # Approx for pi. (3.0, 4.5% error)
    g_frac = (9, 1)           # Approx for g. (9.0, 8% error)
    
    print("\nSelected Titan fractions:")
    print(f"d/h          = {dist_ratio_frac[0]}/{dist_ratio_frac[1]}")
    print(f"4/3          = {four_thirds_frac[0]}/{four_thirds_frac[1]}")
    print(f"g            ≈ {g_frac[0]}/{g_frac[1]}")
    print(f"pi           ≈ {pi_frac[0]}/{pi_frac[1]}")
    print(f"(rho * r^3)  ≈ {mass_factor_frac[0]}/{mass_factor_frac[1]}")
    
    # --- Step 3: Fractional Arithmetic Calculation ---
    print("\n### Step 3: Performing the Calculation ###")
    print("The multiplication order is critical to avoid overflow. The chosen order is:")
    print("F = ( ( ( (d/h) * (4/3) ) * (rho*r^3) ) * pi ) * g\n")
    
    def gcd(a, b):
        while b:
            a, b = b, a % b
        return a

    def titan_mul(frac1, frac2):
        n1, d1 = frac1
        n2, d2 = frac2
        
        # Cross-simplify before multiplication to prevent overflow
        common1 = gcd(n1, d2)
        n1 //= common1
        d2 //= common1
        
        common2 = gcd(n2, d1)
        n2 //= common2
        d1 //= common2
        
        num = n1 * n2
        den = d1 * d2
        
        if num > 31 or den > 31:
            raise ValueError(f"Overflow calculating ({frac1[0]}/{frac1[1]}) * ({frac2[0]}/{frac2[1]}) -> {num}/{den}")
            
        return (num, den)

    # Perform calculation step-by-step
    try:
        # Step A: (d/h) * (4/3)
        res1 = titan_mul(dist_ratio_frac, four_thirds_frac)
        print(f"Step A: ({dist_ratio_frac[0]}/{dist_ratio_frac[1]}) * ({four_thirds_frac[0]}/{four_thirds_frac[1]}) = {res1[0]}/{res1[1]}")
        
        # Step B: res1 * (rho*r^3)
        res2 = titan_mul(res1, mass_factor_frac)
        print(f"Step B: ({res1[0]}/{res1[1]}) * ({mass_factor_frac[0]}/{mass_factor_frac[1]}) = {res2[0]}/{res2[1]}")
        
        # Step C: res2 * pi
        res3 = titan_mul(res2, pi_frac)
        print(f"Step C: ({res2[0]}/{res2[1]}) * ({pi_frac[0]}/{pi_frac[1]}) = {res3[0]}/{res3[1]}")
        
        # Step D: res3 * g
        final_result_frac = titan_mul(res3, g_frac)
        print(f"Step D: ({res3[0]}/{res3[1]}) * ({g_frac[0]}/{g_frac[1]}) = {final_result_frac[0]}/{final_result_frac[1]}")

        # --- Step 4: Final Result and Error Analysis ---
        print("\n" + "-" * 30)
        print("\n### Step 4: Final Result and Error ###")

        print("\nThe final equation with the chosen approximations is:")
        print(f"F = ({dist_ratio_frac[0]}/{dist_ratio_frac[1]}) * ({g_frac[0]}/{g_frac[1]}) * ({four_thirds_frac[0]}/{four_thirds_frac[1]}) * ({pi_frac[0]}/{pi_frac[1]}) * ({mass_factor_frac[0]}/{mass_factor_frac[1]})")

        final_force_val = final_result_frac[0] / final_result_frac[1]
        print(f"\nThe calculated force is F = {final_result_frac[0]}/{final_result_frac[1]} = {final_force_val:.1f} N")

        true_force_val = (d_si / h_si) * g_si * (4.0/3.0) * pi_val * (rho_si * r_si**3)
        abs_error = abs(true_force_val - final_force_val)
        
        print(f"\nFor reference, the true force value is approximately {true_force_val:.5f} N.")
        print(f"The absolute error of the Titan calculation is {abs_error:.5f}.")
        
    except ValueError as e:
        print(f"\nCalculation failed due to overflow: {e}")

if __name__ == '__main__':
    main()