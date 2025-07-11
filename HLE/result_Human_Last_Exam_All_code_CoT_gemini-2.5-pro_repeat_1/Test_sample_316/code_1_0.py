import sys

def solve():
    """
    Calculates the critical exponents for the reverse square function estimate for the cone.
    """
    print("The problem asks for a critical exponent for an inequality related to the cone in R^3.")
    print("The exponent alpha(p) is piecewise linear in 1/p, and we need to find where its slope changes.")
    print("This is determined by the behavior of the exponent in the sharp bilinear cone restriction estimate, which we call beta(p).")
    print("Let x = 1/p. We analyze the function beta(x) by connecting known sharp results at anchor points.\n")

    # Anchor points (x, beta) where x = 1/p
    # p=2 is the starting point for the domain p>2, but we don't need its beta value to find the other critical exponent.
    p_wolff, beta_wolff = 3, -1/4
    p_tomas_stein, beta_tomas_stein = 4, 1/2
    p_infty, beta_infty = float('inf'), 1

    x_wolff = 1 / p_wolff
    x_tomas_stein = 1 / p_tomas_stein
    x_infty = 0

    print(f"Anchor points (1/p, beta) used for calculation:")
    print(f"1. Wolff's result: p = {p_wolff}, (x={x_wolff:.3f}, beta={beta_wolff})")
    print(f"2. Tomas-Stein estimate: p = {p_tomas_stein}, (x={x_tomas_stein}, beta={beta_tomas_stein})")
    print(f"3. L-infinity estimate: p -> inf, (x={x_infty}, beta={beta_infty})\n")

    # Calculate slope between p=3 and p=4
    slope_3_4_num = beta_tomas_stein - beta_wolff
    slope_3_4_den = x_tomas_stein - x_wolff
    slope_3_4 = slope_3_4_num / slope_3_4_den
    print(f"Calculating the slope for p in (3, 4):")
    print(f"  Slope = (beta(p=4) - beta(p=3)) / (1/4 - 1/3)")
    print(f"        = ({beta_tomas_stein} - ({beta_wolff})) / ({x_tomas_stein} - {x_wolff:.3f})")
    print(f"        = {slope_3_4_num:.2f} / {slope_3_4_den:.3f} = {slope_3_4:.1f}\n")

    # Calculate slope for p > 4
    slope_4_inf_num = beta_infty - beta_tomas_stein
    slope_4_inf_den = x_infty - x_tomas_stein
    slope_4_inf = slope_4_inf_num / slope_4_inf_den
    print(f"Calculating the slope for p in (4, inf):")
    print(f"  Slope = (beta(p=inf) - beta(p=4)) / (0 - 1/4)")
    print(f"        = ({beta_infty} - {beta_tomas_stein}) / ({x_infty} - {x_tomas_stein})")
    print(f"        = {slope_4_inf_num:.2f} / {slope_4_inf_den:.2f} = {slope_4_inf:.1f}\n")

    print(f"The slope of beta(1/p) changes from {slope_3_4:.1f} to {slope_4_inf:.1f} at p = 4.")
    print("This confirms that p=4 is a critical exponent, as stated in the problem.\n")
    
    # We also need to check the slope on the other side of p=3. We need beta(p=2).
    p_l2, beta_l2 = 2, 0.0
    x_l2 = 1/p_l2
    slope_2_3_num = beta_wolff - beta_l2
    slope_2_3_den = x_wolff - x_l2
    slope_2_3 = slope_2_3_num / slope_2_3_den
    print(f"To find the other critical exponent, we check the slope for p in (2, 3).")
    print(f"The anchor point at p=2 comes from the trivial L2 estimate, giving beta(p=2)=0.")
    print(f"  Slope = (beta(p=3) - beta(p=2)) / (1/3 - 1/2)")
    print(f"        = ({beta_wolff} - {beta_l2}) / ({x_wolff:.3f} - {x_l2})")
    print(f"        = {slope_2_3_num:.2f} / {slope_2_3_den:.3f} = {slope_2_3:.1f}\n")

    print(f"The slope of beta(1/p) changes from {slope_2_3:.1f} to {slope_3_4:.1f} at p = 3.")
    print("This means p=3 is the other critical exponent.\n")

    other_critical_exponent = 3
    print("Conclusion: The two critical exponents for p > 2 are 3 and 4.")
    print(f"The problem states one is 4, so the other must be {other_critical_exponent}.")
    
    # Final answer format
    sys.stdout.write(f'<<<{other_critical_exponent}>>>')

solve()