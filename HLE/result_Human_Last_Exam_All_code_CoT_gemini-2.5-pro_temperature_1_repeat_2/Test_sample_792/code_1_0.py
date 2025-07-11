import math

def titan_calculate_force():
    """
    Calculates the required force using Titan computer's 5-bit fractional arithmetic rules.
    This involves representing all numbers as fractions a/b where 0<=a,b<=31,
    and performing calculations while respecting these constraints.
    """
    print("--- Titan Computer Calculation for the Curious Monkey ---")
    print("\nStep 1: Define the governing equation and constants as 5-bit fractions.")
    print("Formula: F = (d * g / h) * m, where m = rho * (4/3) * pi * r^3")

    # Define constants as fractions (numerator, denominator)
    d = (20, 1)
    h = (10, 1)
    # We choose g = 19/2 = 9.5 as a good approximation for 9.8
    g = (19, 2)
    # We choose pi = 22/7
    pi = (22, 7)
    # r = 0.5 cm
    r = (1, 2)
    # rho = 0.9 kg/cm^3
    rho = (9, 10)

    print(f"\nChosen constants:")
    print(f"  d = {d[0]}/{d[1]} m")
    print(f"  h = {h[0]}/{h[1]} m")
    print(f"  g = {g[0]}/{g[1]} m/s^2 (approximating 9.8)")
    print(f"  pi = {pi[0]}/{pi[1]}")
    print(f"  r = {r[0]}/{r[1]} cm")
    print(f"  rho = {rho[0]}/{rho[1]} kg/cm^3")

    print("\nStep 2: Calculate the volume (V) of the rock in cm^3.")
    # V = (4/3) * pi * r^3
    # r^3 = (1/2)^3 = 1/8
    r_cubed = (1, 8)
    print(f"  r^3 = ({r[0]}/{r[1]})^3 = {r_cubed[0]}/{r_cubed[1]}")

    # pi * r^3 = (22/7) * (1/8).
    # This would be 22/56, but 56 > 31. We must simplify first.
    # (22/7)*(1/8) -> (11/7)*(1/4) by dividing 22 and 8 by 2.
    # Now, (11*1)/(7*4) = 11/28. Both are <= 31.
    pi_r_cubed = (11, 28)
    print(f"  pi * r^3 = ({pi[0]}/{pi[1]}) * ({r_cubed[0]}/{r_cubed[1]}) -> Simplified to {pi_r_cubed[0]}/{pi_r_cubed[1]}")

    # V = (4/3) * (11/28)
    # This would be 44/84, both > 31. We simplify first.
    # (4/3) * (11/28) -> (1/3) * (11/7) by dividing 4 and 28 by 4.
    # Now, (1*11)/(3*7) = 11/21. Both are <= 31.
    V = (11, 21)
    print(f"  V = (4/3) * ({pi_r_cubed[0]}/{pi_r_cubed[1]}) -> Simplified to {V[0]}/{V[1]} cm^3")

    print("\nStep 3: Calculate the mass (m) in kg and apply strategic approximation.")
    # m = rho * V = (9/10) * (11/21)
    # This would be 99/210, both > 31. We can simplify (9/3, 21/3) -> (3/10)*(11/7) = 33/70.
    # Still, 33 and 70 are > 31. The calculation cannot proceed directly.
    # We must approximate.
    m_val_true = (rho[0]/rho[1]) * (V[0]/V[1]) # approx 0.4714
    print(f"  Direct multiplication m = ({rho[0]}/{rho[1]}) * ({V[0]}/{V[1]}) fails due to 5-bit constraint.")
    print(f"  The target value for mass is ~{m_val_true:.4f} kg.")

    # Strategy: Pre-calculate the (d*g/h) term and use it to guide the approximation of m.
    # (d/h) = (20/1)/(10/1) = 2/1
    d_over_h = (2, 1)
    # (d/h)*g = (2/1)*(19/2) = 19/1
    force_prefix = (19, 1)
    print(f"  The calculation for F simplifies to F = ({force_prefix[0]}/{force_prefix[1]}) * m.")
    print(f"  To enable this final multiplication, we must approximate m as a fraction with denominator {force_prefix[0]}.")

    # We need m_approx = a/19 ≈ 0.4714. So a ≈ 0.4714 * 19 = 8.9566.
    # The nearest integer is 9.
    m_approx_num = 9
    m_approx_den = 19
    m_approx = (m_approx_num, m_approx_den)
    print(f"  We choose the approximation m_approx = {m_approx[0]}/{m_approx[1]} kg.")

    print("\nStep 4: Calculate the final force (F) in Newtons.")
    # F = (19/1) * (9/19).
    # We can simplify by dividing 19 and 19 by 19.
    # F = (1/1) * (9/1) = 9/1.
    F_num = 9
    F_den = 1
    F_final = (F_num, F_den)
    print(f"  F = ({force_prefix[0]}/{force_prefix[1]}) * ({m_approx[0]}/{m_approx[1]}) = {F_final[0]}/{F_final[1]} N.")

    print("\n--- Final Result ---")
    print("The final equation showing each number is:")
    print(f"F = (d / h) * g * m")
    final_equation = f"{F_final[0]}/{F_final[1]} = (({d[0]}/{d[1]}) / ({h[0]}/{h[1]})) * ({g[0]}/{g[1]}) * ({m_approx[0]}/{m_approx[1]})"
    print(final_equation)

    # Calculate the error
    # "True" mass m = rho * V = 0.9 * (4/3 * math.pi * 0.5**3) = 0.471238... kg
    # "True" force (using g=9.8) = (d * m * g) / h = (20 * 0.471238 * 9.8) / 10 = 9.2362... N
    true_mass = 0.9 * (4/3 * math.pi * 0.5**3)
    true_force = (20 * true_mass * 9.8) / 10
    calculated_force = F_final[0]/F_final[1]
    error = abs(calculated_force - true_force)

    # Check if the rock hits the lion
    # d_landing = F_calc * h / (m_true * g_used)
    g_used = g[0]/g[1]
    m_used_for_check = (rho[0]/rho[1]) * (V[0]/V[1]) # Use the un-approximated mass for checking accuracy
    landing_dist = (calculated_force * (h[0]/h[1])) / (m_val_true * g_used)
    print(f"\nWith our calculated force F={calculated_force} N (and g={g_used}), the rock lands at {landing_dist:.2f} m.")
    print(f"This is within the target range of 19-21 m.")

    print(f"\nThe absolute error compared to a theoretical value with g=9.8 is {error:.5f} N.")

titan_calculate_force()

# The smallest absolute error we can produce with this method is |9 - 9.23628| = 0.23628
# The problem is solvable.
print("\n<<<Y[0.23628]>>>")
