import math

def solve_shock_profile():
    """
    Calculates and prints the analytical solution for a shock wave density profile.

    This function follows the theoretical derivation for a Mach sqrt(3) shock in a
    monatomic ideal gas with a Prandtl number of 3/4. It presents the density
    profile rho(xi) / rho_0, where xi = x / L is the non-dimensional position.
    """

    # 1. Define the given physical constants
    M_sq = 3.0
    M = math.sqrt(M_sq)
    gamma = 5.0 / 3.0  # For a monatomic ideal gas

    # 2. Use Rankine-Hugoniot relations to find the density/velocity ratio across the shock.
    # rho_1 / rho_0 = u_0 / u_1 = ((gamma + 1) * M^2) / ((gamma - 1) * M^2 + 2)
    density_ratio = ((gamma + 1) * M_sq) / ((gamma - 1) * M_sq + 2)
    velocity_ratio = 1.0 / density_ratio # u_1 / u_0

    # The velocity profile u(x) has the form of a hyperbolic tangent:
    # u(x)/u_0 = ( (1+u_1/u_0)/2 - (1-u_1/u_0)/2 * tanh(C * xi) )
    # From this, the density profile rho(x)/rho_0 = u_0/u(x) can be derived.
    # For u_1/u_0 = 1/2, the expression simplifies to:
    # rho(xi)/rho_0 = 4 / (3 - tanh(C * xi))
    
    # 3. Determine the coefficient C in the tanh argument.
    # The argument is x / (2 * delta), where delta is the characteristic shock thickness.
    # delta is derived from the full Navier-Stokes equations for Pr=3/4.
    # The dimensionless argument C * xi is given by (gamma * M * c_0 / 4) * xi.
    # Assuming a non-dimensional system where the ambient sound speed c_0 = 1.
    C = (gamma * M) / 4.0

    # 4. Print the final analytical solution step-by-step.
    print("The analytical solution for the density profile ρ(ξ)/ρ₀ is derived as follows:")
    print("where ξ = x/L is the position in units of the ambient conductive length scale L.")
    print("-" * 70)
    
    print(f"Given parameters:")
    print(f"  Shock Mach number, M = sqrt({M_sq}) ≈ {M:.4f}")
    print(f"  Specific heat ratio for monatomic gas, γ = {gamma:.4f}")
    print("\nStep 1: Determine the post-shock to pre-shock density ratio (ρ₁/ρ₀).")
    print(f"  ρ₁/ρ₀ = ( (γ+1)M² ) / ( (γ-1)M² + 2 )")
    print(f"  ρ₁/ρ₀ = ( ({gamma:.2f}+1)*{M_sq:.1f} ) / ( ({gamma:.2f}-1)*{M_sq:.1f} + 2 ) = {density_ratio:.1f}")
    
    print("\nStep 2: State the general form of the Becker-Taylor shock profile solution.")
    print("  For Pr = 3/4, the velocity profile is a tanh function, leading to a density profile:")
    print("  ρ(ξ)/ρ₀ = 4 / (3 - tanh(C * ξ))")
    
    print("\nStep 3: Calculate the dimensionless coefficient C for the given parameters.")
    print("  The coefficient C in the argument of the tanh function is C = (γ * M) / 4.")
    print(f"  C = ({gamma:.4f} * {M:.4f}) / 4.0 ≈ {C:.4f}")
    
    print("-" * 70)
    print("Final Analytical Solution:")
    print("The density profile is given by the equation:")
    print("\n  ρ(ξ)     4")
    print("------ = ------------------------")
    print(f"  ρ₀      3 - tanh({C:.4f} * ξ)")
    print("-" * 70)

solve_shock_profile()