import numpy as np
from scipy.special import airy

def solve():
    """
    Calculates the radius of the spherical balloon y(t) at t=pi/4.
    
    The solution is based on the hypothesis that y(t) is an Airy function,
    as suggested by the initial conditions. The hypothesized solution that
    matches the initial conditions y(0) and y'(0)=0 is:
    y(t) = (1/128) * sec(t) * Bi(-t**4)
    """
    
    # Value of t
    t = np.pi / 4
    
    # Calculate sec(t)
    sec_t = 1 / np.cos(t)
    
    # The argument for the Airy function Bi is x = -t**4
    x = -(t**4)
    
    # Use scipy.special.airy to get the value of Bi(x).
    # airy(x) returns a tuple (Ai, Aip, Bi, Bip)
    ai, aip, bi, bip = airy(x)
    
    # Calculate the final radius y(pi/4)
    y_pi_over_4 = (1/128) * sec_t * bi
    
    # Print the calculation steps
    print(f"The radius is to be found at t = pi/4.")
    print(f"The hypothesized solution is y(t) = (1/128) * sec(t) * Bi(-t^4).")
    print(f"We evaluate this at t = pi/4 = {t:.6f}")
    print(f"The equation to compute is: y(pi/4) = (1/128) * sec(pi/4) * Bi(-(pi/4)^4)")
    print(f"1. sec(pi/4) = {sec_t:.6f}")
    print(f"2. Argument for Bi: x = -(pi/4)^4 = {x:.6f}")
    print(f"3. Bi(x) = Bi({x:.6f}) = {bi:.6f}")
    print(f"4. Final calculation: y(pi/4) = (1/128) * {sec_t:.6f} * {bi:.6f}")
    print(f"Result: y(pi/4) = {y_pi_over_4:.6f}")

solve()