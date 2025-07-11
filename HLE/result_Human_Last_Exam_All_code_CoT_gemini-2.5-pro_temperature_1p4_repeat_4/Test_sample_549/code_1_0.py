import sympy
from sympy import integrate, Symbol, pi, sqrt, Rational

def solve_conductivity_correction():
    """
    This script evaluates the quantum correction to conductivity for an electron
    in a 3D bulk semiconductor based on the principle of weak localization.
    """
    # Introduction and plan
    print("This script evaluates the quantum correction to conductivity (Δσ) for an electron in a 3D bulk semiconductor.")
    print("Based on the theory of weak localization, this correction is proportional to the probability of the electron returning to its origin.")
    print("We will calculate this return probability, W, by integrating the return probability density P(0, t) over the relevant time interval.\n")

    # Define symbolic variables for the calculation
    # D: Diffusion coefficient
    # t: time
    # τ: Mean free time (lower integration limit)
    # τ_φ: Phase coherence time (upper integration limit)
    D = Symbol('D', positive=True, real=True)
    t = Symbol('t', positive=True, real=True)
    tau = Symbol('τ', positive=True, real=True)
    tau_phi = Symbol('τ_φ', positive=True, real=True)

    # Step 1: Define the return probability density in 3D
    print("Step 1: The probability density for a diffusing particle to be at its origin at time t in 3 dimensions is given by:")
    # P(0, t) = (4 * pi * D * t)^(-3/2)
    P_0_t = (4 * pi * D * t)**(-Rational(3, 2))
    print("P(0, t) = (4 ⋅ π ⋅ D ⋅ t)^(-3/2)")
    print("-" * 50)

    # Step 2: Integrate P(0, t) from τ to τ_φ to find the total return probability W
    print(f"Step 2: We integrate P(0, t) from the mean free time τ to the phase coherence time τ_φ.")
    print("W = ∫ P(0, t) dt  (from t=τ to t=τ_φ)\n")
    
    # Perform the integration
    W = integrate(P_0_t, (t, tau, tau_phi))
    
    # Simplify the expression for better readability
    # The integral of t^(-3/2) is -2*t^(-1/2)
    # The (4*pi*D)^(-3/2) factor is a constant
    # Result is (4*pi*D)^(-3/2) * [-2*t^(-1/2)] from τ to τ_φ
    # = (4*pi*D)^(-3/2) * 2 * (1/sqrt(τ) - 1/sqrt(τ_φ))
    # = (1/(8*pi^(3/2)*D^(3/2))) * 2 * (1/sqrt(τ) - 1/sqrt(τ_φ))
    # = (1/(4*pi^(3/2)*D^(3/2))) * (1/sqrt(τ) - 1/sqrt(τ_φ))
    
    print("The result for the total return probability W is:")
    print("W = (1/(4 ⋅ π^(3/2) ⋅ D^(3/2))) ⋅ (1/√τ - 1/√τ_φ)")
    print("-" * 50)

    # Step 3: Relate W to the quantum correction Δσ
    print("Step 3: The quantum correction to conductivity, Δσ, is proportional to this probability W.")
    print("Therefore, the final expression for the quantum correction is:\n")

    # Display the final equation clearly
    print("Δσ = C_σ ⋅ W")
    print("where C_σ is a constant of proportionality. Substituting W, we get:")
    
    # Print the final equation with all its numbers and symbols
    # Δσ = C_σ * [ (1/sqrt(τ) - 1/sqrt(τ_φ)) / (4 * π^(3/2) * D^(3/2)) ]
    print("\n   Final Equation:")
    print("   -----------------------------------------------------------")
    print("          C_σ      /    1          1      \\")
    print("   Δσ = --------- ⋅ | -------  -  ------- |")
    print("        4π³ D³   \\  √τ         √τ_φ   /")
    # A cleaner text-based representation
    print("\nOr, in a single line:")
    print("Δσ = C_σ * [ (1/sqrt(τ) - 1/sqrt(τ_φ)) / (4 * π^(3/2) * D^(3/2)) ]")
    print("   -----------------------------------------------------------")


if __name__ == '__main__':
    solve_conductivity_correction()