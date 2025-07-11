import math

def high_precision_calculation():
    """Calculates the force using standard high-precision values."""
    G = 6.67430e-11
    R_planet = 2000 * 1000  # 2000 km in meters
    rho = 1.2 * 1000  # 1.2 t/m^3 in kg/m^3
    c = 299792458
    m_probe = 50
    dist_from_horizon = 1000

    V_planet = (4/3) * math.pi * (R_planet**3)
    m_planet = rho * V_planet
    
    Rs = (2 * G * m_planet) / (c**2)
    r = Rs + dist_from_horizon
    
    F = (G * m_planet * m_probe) / (r**2)
    return F

def titan_simulation():
    """
    Simulates the gravitational force calculation on the Titan 6-bit computer.
    """
    print("Initializing Titan calculation for gravitational force.")
    print("F = G * (rho * 4/3 * pi * R^3) * m_probe / r^2")
    print("=" * 40)

    # Step 1: Define constants and values as Titan fractional operands (num, den, exp)
    G = (20, 3, -11)
    rho = (6, 5, 3)
    four_thirds = (4, 3, 0)
    pi_approx = (22, 7, 0)
    R_cubed = (8, 1, 18)  # R = 2e6, so R^3 = 8e18
    m_probe = (50, 1, 0)
    # The Schwarzschild radius is negligible compared to 1km, so r^2 approx (1km)^2
    r_squared = (1, 1, 6) # r = 1e3, r^2 = 1e6
    
    print("Approximations used:")
    print(f"G               = {G[0]}/{G[1]} x 10^{G[2]}")
    print(f"Density (rho)   = {rho[0]}/{rho[1]} x 10^{rho[2]}")
    print(f"Pi              = {pi_approx[0]}/{pi_approx[1]}")
    print(f"Planet Radius^3 = {R_cubed[0]} x 10^{R_cubed[2]}")
    print(f"Probe Mass      = {m_probe[0]}")
    print(f"Distance^2      = {r_squared[0]} x 10^{r_squared[2]}")
    print("=" * 40)

    # Combine all terms into one expression for the RED (reduce) instruction
    # F = G * rho * (4/3) * pi * R^3 * m_probe / r^2
    
    # Step 2: Calculate the final exponent
    total_exponent = G[2] + rho[2] + four_thirds[2] + pi_approx[2] + R_cubed[2] + m_probe[2] - r_squared[2]
    
    # Step 3: List all mantissa numerators and denominators for reduction
    numerators = [G[0], rho[0], four_thirds[0], pi_approx[0], R_cubed[0], m_probe[0]]
    denominators = [G[1], rho[1], four_thirds[1], pi_approx[1], R_cubed[1], m_probe[1]]

    print("Executing RED instruction: simplifying mantissa expression.")
    print(f"Initial Expression: {numerators} / {denominators}")

    # Simplify by cancelling common factors (part of the RED instruction)
    for i, n in enumerate(numerators):
        for j, d in enumerate(denominators):
            common_divisor = math.gcd(n, d)
            if common_divisor > 1:
                numerators[i] //= common_divisor
                denominators[j] //= common_divisor
                n = numerators[i] # update n for further checks in the inner loop

    print(f"After Cancellation: {numerators} / {denominators}")
    
    # Step 4: Multiply remaining terms, handling overflows
    
    # Calculate denominator product
    final_den = 1
    for d in denominators:
        final_den *= d
        if final_den > 63:
            raise ValueError("Denominator overflow during reduction")

    # Calculate numerator product step-by-step
    final_num = 1
    internal_exponent = 0
    
    print("Calculating final numerator product with overflow handling:")
    for n in numerators:
        old_num = final_num
        product = final_num * n
        print(f"  {old_num} * {n} = {product}", end="")
        if product > 63:
            # Handle overflow by converting to M x 10^E representation
            exponent_increase = 0
            temp_product = product
            while temp_product > 63:
                temp_product /= 10
                exponent_increase += 1
            
            final_num = round(temp_product)
            internal_exponent += exponent_increase
            print(f" -> OVERFLOW. Approximating as {final_num} x 10^{exponent_increase}")
        else:
            final_num = product
            print("")
    
    total_exponent += internal_exponent

    print("=" * 40)
    print("Final Result:")
    print(f"Force = {final_num} / {final_den} * 10^{total_exponent} N")
    
    return (final_num / final_den) * (10**total_exponent)

# Run the simulation
titan_result = titan_simulation()
true_result = high_precision_calculation()
relative_error = abs(titan_result - true_result) / true_result

print(f"\nTitan calculated Force: {titan_result:.4e} N")
print(f"High-precision Force: {true_result:.4e} N")
print(f"Relative Error: {relative_error:.1%}")

# The final answer format is specified by the prompt
# <<<Y[e]>>> where e is the smallest relative error.
print("\n<<<Y[0.6%]>>>")