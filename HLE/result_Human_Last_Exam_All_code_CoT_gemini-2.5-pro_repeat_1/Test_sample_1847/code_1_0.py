import math

def display_solution():
    """
    Prints the analytical solution for the density profile of a Mach sqrt(3) shock wave.
    """
    
    # The derived analytical solution for the density profile rho(x).
    # It is expressed in terms of normalized density (rho/rho_0) and normalized position (x/L).
    
    print("The analytical solution for the density profile ρ(x) is given by the following equation:")
    print("")
    
    # Printing the equation with Unicode characters for better formatting
    print("  ρ(x)          _______________________________________________________")
    print("  ────  = 1 +  ╱         4x/3L                4x/3L                   ")
    print("   ρ₀         ╱   exp( ─────── ) * ( 8 + exp( ─────── ) )            ")
    print("             ╱   ─────────────────────────────────────────────────────")
    print("            ╱                             4x/3L                       ")
    print("          ╲╱                        8 + exp( ─────── )                ")
    print("")
    print("where:")
    print("  ρ(x) : Gas density at position x.")
    print("  ρ₀   : Ambient (pre-shock) density.")
    print("  x    : Position coordinate in the direction of shock propagation.")
    print("  L    : Ambient conductive length scale, defined as L = κ / (ρ₀ * M * c₀ * Cv).")
    print("         - κ: Thermal conductivity")
    print("         - M: Mach number (√3 for this problem)")
    print("         - c₀: Ambient speed of sound")
    print("         - Cv: Specific heat at constant volume")
    print("\nThe numerical constants in the equation are 8, 1, 4, and 3.")


if __name__ == "__main__":
    display_solution()