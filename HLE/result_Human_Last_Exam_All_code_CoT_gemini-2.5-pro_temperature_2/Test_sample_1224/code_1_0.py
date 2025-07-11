import math

def calculate_and_print_landing_time():
    """
    Calculates and explains the landing time calculation on the Titan 4-bit computer.
    """
    print("### Titan Landing Time Calculation for Pandora ###")
    print("\nStep 1: State formulas and approximations.")
    print("The landing time 't' is calculated using the formula t = sqrt(2 * h / g).")
    print("With a simplified spherical model for Pandora, g = (4/3) * pi * G * R * d.")
    print("This leads to t^2 = (3 * h) / (2 * pi * G * R * d).")
    
    print("\nTo stay within Titan's 4-bit constraints, we use these approximations:")
    
    # Values for the equation
    val_3 = 3
    val_h = 5000
    val_2 = 2
    val_pi = 3
    val_G_str = "5e-11"
    val_G_num = 5e-11
    val_R = 2e6
    val_d = 300
    
    # Printing approximations as they are used in the Titan architecture
    print(f"- h = {val_h} m -> 5/1 * 10^3")
    print(f"- pi approx = {val_pi}/1")
    print(f"- G approx = {val_G_str} N m^2/kg^2 -> 5/1 * 10^-11")
    print(f"- R = {int(val_R)} m -> 2/1 * 10^6")
    print(f"- d = {val_d} kg/m^3 -> 3/1 * 10^2")

    print("\nStep 2: Present the full equation for t^2 with approximated values.")
    print("t^2 = ({} * {}) / ({} * {} * {} * {} * {})".format(
        val_3, val_h, val_2, val_pi, val_G_str, int(val_R), val_d
    ))

    # Perform the calculation
    t_squared_num = val_3 * val_h
    t_squared_den = val_2 * val_pi * val_G_num * val_R * val_d
    t_squared_val = t_squared_num / t_squared_den

    print("\nStep 3: Show simplification of the t^2 expression.")
    print(f"This expression simplifies to t^2 = {int(t_squared_num * 10000)} / {int(t_squared_den * 10000)}.")
    # On Titan, this is computed sequentially to avoid overflow. The result is 1/12 * 10^6
    print(f"On Titan, this results in a fraction and exponent: 1/12 * 10^6.")
    
    print("\nStep 4: Calculate the square root using pre-computed approximation.")
    print(f"t = sqrt(1/12 * 10^6) = sqrt(1/12) * 10^3.")
    # Sqrt(1/12) is approx 0.2887. Best 4-bit fraction is 2/7 = 0.2857.
    sqrt_approx_num = 2
    sqrt_approx_den = 7
    print(f"Rule 3a allows using pre-computed approximations for functions like sqrt.")
    print(f"We use the 4-bit friendly fraction {sqrt_approx_num}/{sqrt_approx_den} to approximate sqrt(1/12).")
    
    print(f"\nFinal Titan Calculation: t = {sqrt_approx_num}/{sqrt_approx_den} * 10^3")
    
    # Final value from Titan
    titan_t = (sqrt_approx_num / sqrt_approx_den) * 1000

    print("\n-----------------------------------------")
    print(f"Final landing time calculated by Titan: {titan_t:.2f} seconds.")
    print("-----------------------------------------")

calculate_and_print_landing_time()

# Benchmark calculation for error
# M = (4/3)*pi*[rc^3*(dc-ds) + a^2*c*ds]
M_true = (4/3)*math.pi * ( (1e5)**3 * (1200-300) + (2e6)**2 * (1.985e6) * 300 )
g_true = 6.6743e-11 * M_true / (2e6)**2
t_true = math.sqrt(2 * 5000 / g_true)
error = abs(285.71 - t_true)
# print(f"Benchmark t={t_true}, titan t={285.71}, error={error}")
# <<<Y[40.71]>>>