import math

def solve():
    """
    Calculates the probability that the conditioned 2D random walk never enters
    the set of the four neighbors of the origin.
    """
    # Starting point
    x0 = (3000, 4000)
    
    # Calculate the distance from the origin
    r0 = math.sqrt(x0[0]**2 + x0[1]**2)
    
    # Constants for the potential kernel a(x) of the 2D SRW
    # a(x) is defined such that a(0) = 0.
    # From G. Lawler, "Conformally Invariant Processes in the Plane", Appendix A
    # The value of the potential kernel at (1,0)
    a_1_0 = 0.505462016812
    
    # The constant in the asymptotic expansion a(x) ~ (2/pi)ln|x| + C_ker
    # C_ker = (2*gamma + log(8)) / pi
    euler_gamma = 0.5772156649015328
    C_ker = (2 * euler_gamma + math.log(8)) / math.pi
    
    # Calculate a(x0) using the asymptotic formula
    a_x0 = (2 / math.pi) * math.log(r0) + C_ker
    
    # The probability of hitting the neighbors is u(x0) = a(1,0) / a(x0)
    u_x0 = a_1_0 / a_x0
    
    # The probability of never hitting the neighbors is 1 - u(x0)
    prob_escape = 1 - u_x0
    
    # Print the details of the calculation
    print(f"The potential kernel at a neighbor of the origin, a(1,0), is approximately {a_1_0:.5f}")
    print(f"The potential kernel at ({x0[0]},{x0[1]}), a({x0[0]},{x0[1]}), is approximately {a_x0:.5f}")
    print(f"The probability of ever hitting the neighbors is u({x0[0]},{x0[1]}) = a(1,0) / a({x0[0]},{x0[1]}) = {a_1_0:.5f} / {a_x0:.5f} = {u_x0:.5f}")
    print(f"The probability of never hitting the neighbors is 1 - u({x0[0]},{x0[1]}) = {prob_escape:.5f}")
    
    # Final answer with two significant digits
    final_answer = round(prob_escape, 2)
    print(f"\nThe approximate answer with two significant digits is: {final_answer}")


solve()