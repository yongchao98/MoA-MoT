import math

def solve_radiation_ratio():
    """
    Calculates the maximum achievable ratio of bidirectional conical power
    to line intensity for the described rotating charged sphere.

    The derivation shows the optimal radiation pattern is I(θ) ∝ 1 + cos²(θ).
    The ratio is maximized by taking the line intensity at the minimum (θ=π/2)
    and the cone power centered on the maximum (θ=0 and θ=π).

    The final ratio is given by the formula:
    Ratio = 2 * (Integral of I(θ) over one cone) / (I_min)
          = 2 * [∫(φ=0 to 2π) dφ ∫(θ=0 to π/4) (1+cos²θ)sinθ dθ] / [1+cos²(π/2)]
          = 4π * [ -cosθ - (cos³θ)/3 ] from θ=0 to π/4
          = 4π * ( ( -1/√2 - 1/(6√2) ) - ( -1 - 1/3 ) )
          = 4π * ( 4/3 - 7/(6√2) )
          = (π/3) * (16 - 7√2)
    """

    # Numbers in the final symbolic equation
    pi_val = math.pi
    sqrt2_val = math.sqrt(2)
    c1 = 16
    c2 = 7
    c3 = 3

    # Calculate the final result
    result = (pi_val / c3) * (c1 - c2 * sqrt2_val)

    # Print the explanation and calculation
    print("The final equation for the maximum achievable ratio is:")
    print(f"Ratio = (π / {c3}) * ({c1} - {c2} * √2)")
    print("\nWhere the values are:")
    print(f"π ≈ {pi_val}")
    print(f"√2 ≈ {sqrt2_val}")
    
    print(f"\nSubstituting the values:")
    print(f"Ratio = ({pi_val:.4f} / {c3}) * ({c1} - {c2} * {sqrt2_val:.4f})")
    print(f"Ratio = {pi_val / c3:.4f} * ({c1} - {c2 * sqrt2_val:.4f})")
    print(f"Ratio = {pi_val / c3:.4f} * {c1 - c2 * sqrt2_val:.4f}")
    
    print(f"\nThe maximum ratio is: {result}")

solve_radiation_ratio()
<<<6.388339031582239>>>