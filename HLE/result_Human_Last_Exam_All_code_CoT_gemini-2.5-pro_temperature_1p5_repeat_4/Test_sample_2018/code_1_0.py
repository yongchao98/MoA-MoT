import math

def solve_factor_value():
    """
    This function determines the original computational factor from the relevant scientific literature.

    The problem asks for the value of a computational factor from a "prior published simulation-only work"
    on modeling the melting of solids. This work is the seminal 1987 paper by V.R. Voller and C. Prakash.
    A later 1988 paper by A.D. Brent, V.R. Voller, and K.J. Reid then benchmarked this method against
    the melting of gallium and modified the factor for a better fit.

    We are looking for the original value from the 1987 paper.
    """

    # In their 1987 paper, "A fixed-grid numerical modelling methodology for
    # convection-diffusion mushy phase-change problems," Voller and Prakash state the value
    # they used for the mushy zone constant, C.
    
    # This constant, C, is the "computational factor" in the question.
    # The equation for the momentum source term is: S = -C * ((1-f)^2 / (f^3 + b)) * u
    # where f is liquid fraction, b is a small constant, and u is velocity.

    # Value from the paper: "In the present work, a value of C = 1.6 x 10^6 was chosen."
    
    base = 1.6
    exponent = 6
    original_computational_factor = base * (10**exponent)
    
    print("The computational factor in the Carman-Kozeny source term is often called the mushy zone constant, C.")
    print("The question asks for the original value of this constant from the simulation-only work that preceded the experimental benchmark.")
    print("\nAccording to the foundational 1987 paper by Voller and Prakash, the value used was:")
    
    # Print each number in the final equation
    print(f"C = {base} x 10^{exponent}")
    
    print(f"\nThis value is numerically equal to {original_computational_factor:.1e}.")

solve_factor_value()