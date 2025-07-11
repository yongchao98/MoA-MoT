import sys

def solve_poynting_vector():
    """
    This function calculates and prints the symbolic expression for the Poynting vector
    both inside and outside the described cylindrical rod.
    
    The final expression is broken down into its radial (r_hat) and axial (z_hat) components.
    """
    
    # --- Symbolic variable definitions for clarity in the final printout ---
    # ρ: uniform volume charge density
    # R: radius of the rod
    # E: magnitude of the uniform external electric field
    # v: speed of the rod
    # r: radial distance from the center of the rod
    # ε₀: permittivity of free space
    
    # --- Results of the derivation ---
    # The Poynting vector S is calculated as S = (1/μ₀) * (E_total x B)
    # The final result does not depend on μ₀ because B is proportional to μ₀.
    
    # --- Inside the rod (r < R) ---
    # E_total = (ρ*r / (2*ε₀)) r_hat + E z_hat
    # B = (μ₀*ρ*v*r / 2) φ_hat
    # S_inside = (1/μ₀) * (E_total x B)
    
    # The r_hat component of S comes from (E z_hat) x (B φ_hat)
    S_inside_r_comp_numerator = "E * ρ * v * r"
    S_inside_r_comp_denominator = "2"
    
    # The z_hat component of S comes from (E_rod r_hat) x (B φ_hat)
    S_inside_z_comp_numerator = "ρ**2 * v * r**2"
    S_inside_z_comp_denominator = "4 * ε₀"

    # --- Outside the rod (r > R) ---
    # E_total = (ρ*R**2 / (2*ε₀*r)) r_hat + E z_hat
    # B = (μ₀*ρ*v*R**2 / (2*r)) φ_hat
    # S_outside = (1/μ₀) * (E_total x B)
    
    # The r_hat component of S
    S_outside_r_comp_numerator = "E * ρ * v * R**2"
    S_outside_r_comp_denominator = "2 * r"
    
    # The z_hat component of S
    S_outside_z_comp_numerator = "ρ**2 * v * R**4"
    S_outside_z_comp_denominator = "4 * ε₀ * r**2"
    
    # --- Print the final results ---
    print("The Poynting vector S has two components: a radial one (r_hat) and an axial one (z_hat).")
    print("The expression depends on whether you are inside or outside the rod.")
    print("-" * 70)
    
    print("Inside the rod (for r < R):")
    # Note the minus sign for the radial component from the cross product z_hat x φ_hat = -r_hat
    print(f"S_inside = - ( {S_inside_r_comp_numerator} / {S_inside_r_comp_denominator} ) r_hat   +   ( {S_inside_z_comp_numerator} / {S_inside_z_comp_denominator} ) z_hat")
    
    print("\n" + "-" * 70 + "\n")
    
    print("Outside the rod (for r > R):")
    # Note the minus sign for the radial component
    print(f"S_outside = - ( {S_outside_r_comp_numerator} / {S_outside_r_comp_denominator} ) r_hat   +   ( {S_outside_z_comp_numerator} / {S_outside_z_comp_denominator} ) z_hat")

if __name__ == '__main__':
    solve_poynting_vector()