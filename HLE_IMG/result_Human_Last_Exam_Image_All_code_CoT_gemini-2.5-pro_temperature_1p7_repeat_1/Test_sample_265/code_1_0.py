import math

def solve_weld_stress():
    """
    Calculates the required weld leg size for a C-shaped weld under eccentric loading.
    """
    # 1. Given Parameters
    P = 90 * 1000  # Force in N
    b = 50.0  # Length of horizontal welds in mm
    d = 100.0 # Length of vertical weld in mm
    e_prime = 150.0 # Eccentricity from the end of the cantilever to the weld
    tau_per = 250.0 # Permissible shear stress in N/mm^2

    print("--- 1. Given Data ---")
    print(f"Force (P): {P/1000} kN")
    print(f"Horizontal weld length (b): {b} mm")
    print(f"Vertical weld length (d): {d} mm")
    print(f"Eccentric distance (from beam end): {e_prime} mm")
    print(f"Permissible shear stress (τ_per): {tau_per} N/mm^2\n")

    # 2. Weld Group Properties
    # Using a coordinate system with the origin at the center of the vertical weld.
    # Total length of the weld
    L = 2 * b + d

    # Centroid calculation
    # By symmetry, y_bar = 0.
    # x_bar = (A1*x1 + A2*x2 + A3*x3) / (A1+A2+A3), where A is length here
    x_bar = (b * (b / 2) + b * (b / 2) + d * 0) / L
    
    # Eccentricity (distance from centroid to force)
    e = e_prime + b - x_bar

    # Twisting moment (Torque)
    T = P * e

    print("--- 2. Weld Group Properties ---")
    print(f"Total weld length (L): 2 * {b} + {d} = {L} mm")
    print(f"Centroid location (x_bar from vertical weld): ({b} * {b/2} + {b} * {b/2}) / {L} = {x_bar:.2f} mm")
    print(f"Eccentricity (e): {e_prime} + {b} - {x_bar:.2f} = {e:.2f} mm")
    print(f"Twisting Moment (T): {P/1000} * {e:.2f} = {T/1e6:.4f} kN-m\n")

    # 3. Polar Moment of Inertia (per unit throat thickness, Ju)
    # I_xx = Sum(I_cxx + A*dy^2) for each weld segment
    # For horizontal welds (W1, W2): I_cxx is negligible. A=b, dy=d/2
    # For vertical weld (W3): I_cxx = d^3/12. A=d, dy=0
    I_xx_u = 2 * (b * (d / 2)**2) + (d**3 / 12)
    
    # I_yy = Sum(I_cyy + A*dx^2) for each weld segment
    # For horizontal welds (W1, W2): I_cyy=b^3/12. A=b, dx = b/2 - x_bar
    # For vertical weld (W3): I_cyy is negligible. A=d, dx = -x_bar
    I_yy_u = 2 * ((b**3 / 12) + b * (b / 2 - x_bar)**2) + d * (-x_bar)**2
    
    # Polar moment of inertia per unit throat thickness
    J_u = I_xx_u + I_yy_u

    print("--- 3. Moment of Inertia (per unit throat thickness) ---")
    print(f"I_xx_u = 2*({b}*({d/2})^2) + {d}^3/12 = {I_xx_u:.2f} mm^3")
    print(f"I_yy_u = 2*({b}^3/12 + {b}*({b/2} - {x_bar:.2f})^2) + {d}*(-{x_bar:.2f})^2 = {I_yy_u:.2f} mm^3")
    print(f"Polar Moment of Inertia (J_u): {I_xx_u:.2f} + {I_yy_u:.2f} = {J_u:.2f} mm^3\n")

    # 4. Stress Calculation at Critical Point
    # The critical point is the top-right corner, farthest from the centroid.
    # Coordinates of critical point relative to centroid G:
    x_c = b - x_bar
    y_c = d / 2

    # Primary shear stress per unit throat thickness (acts downwards)
    tau_1_u = P / L
    
    # Secondary shear stress components per unit throat thickness
    tau_2x_u = - (T * y_c) / J_u
    tau_2y_u = (T * x_c) / J_u
    
    # Total stress components
    tau_x_total_u = tau_2x_u
    tau_y_total_u = -tau_1_u + tau_2y_u

    # Resultant shear stress per unit throat thickness
    tau_res_u = math.sqrt(tau_x_total_u**2 + tau_y_total_u**2)
    
    print("--- 4. Stresses at Critical Point (per unit throat) ---")
    print(f"Coordinates of critical point relative to centroid: (x_c={x_c:.2f}, y_c={y_c:.2f})")
    print(f"Primary shear stress (τ'_1): {P} / {L} = {tau_1_u:.2f} N/mm")
    print(f"Secondary shear stress components (τ'_2x, τ'_2y):")
    print(f"  τ'_2x = -({T:.0f} * {y_c:.2f}) / {J_u:.2f} = {tau_2x_u:.2f} N/mm")
    print(f"  τ'_2y = ({T:.0f} * {x_c:.2f}) / {J_u:.2f} = {tau_2y_u:.2f} N/mm")
    print(f"Resultant stress (τ'_res): sqrt({tau_x_total_u:.2f}^2 + ({-tau_1_u:.2f} + {tau_2y_u:.2f})^2) = {tau_res_u:.2f} N/mm\n")

    # 5. Determine Weld Size
    # Required throat thickness: t = tau_res_u / tau_per
    t_req = tau_res_u / tau_per
    
    # Required weld leg size (h): h = t / 0.707
    h_req = t_req / 0.707

    print("--- 5. Final Weld Size Calculation ---")
    print("The maximum resultant shear stress is τ_max = τ'_res / t.")
    print("Setting τ_max = τ_permissible, we find the required throat thickness 't'.")
    print(f"t_req = τ'_res / τ_per = {tau_res_u:.2f} / {tau_per:.2f} = {t_req:.3f} mm")
    print("\nThe size of the weld (leg, h) is related to the throat thickness (t) by t = 0.707 * h.")
    print(f"h_req = t_req / 0.707 = {t_req:.3f} / 0.707 = {h_req:.2f} mm")
    print(f"\nThe required size of the weld (leg) is {h_req:.2f} mm.")

    return h_req

if __name__ == '__main__':
    result = solve_weld_stress()
    # The final answer in the required format
    # print(f'<<<{result:.2f}>>>')