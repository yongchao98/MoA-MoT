import numpy as np
from scipy.optimize import root

def get_constant_c():
    """
    This function calculates the smallest possible constant C.
    
    The constant C is found by maximizing a ratio which depends on a parameter alpha.
    The optimal alpha is found by solving a transcendental equation numerically.
    """
    
    # Define the transcendental equation tan(x) = x + pi/2
    def func(x):
        return np.tan(x) - x - np.pi/2

    # Numerically find the root x0 of the equation in the interval (0, pi).
    # A good initial guess is important for the solver. From graphical analysis, 
    # the root is in (pi/2, pi), but the equation can be rewritten for x in (0, pi/2).
    # Our x = 2*alpha, with alpha in [0, pi/2], so x is in [0, pi]. The root is in (0, pi/2).
    # Re-checking the plots, the solution is in (pi/2, 3*pi/2), but my function R(alpha) 
    # considers x = 2alpha where alpha is in [0, pi/2], so x must be in [0, pi].
    # There is indeed a unique solution for x in (0, pi/2). Wait, my derivative analysis leads
    # to tan(x)=x+pi/2, which for x in (0,pi) has a unique solution in (pi/2, pi) if my memory 
    # serves right. Oh wait, my very first graph plot of tan(x) and x+pi/2 had a mistake, I drew x+pi/2 as linear line starting at 0, not at pi/2 on y-axis.
    # In (0, pi/2), tan(x) starts at 0, x+pi/2 starts at pi/2. Both have slope 1 at x=0. Since tan(x) is convex, tan(x) grows faster. So they must cross.
    # So a solution in (0, pi/2) is NOT guaranteed. Let's re-plot.
    #
    # Final check on logic: derivative of R(alpha) is zero when tan(2*alpha) = (pi + 4*alpha)/2.
    # Let x = 2*alpha. tan(x) = pi/2 + x.
    # On (0, pi/2): tan(0)=0, 0+pi/2 = pi/2. tan(x) > x, but maybe not > x+pi/2.
    # g(x) = tan(x) - x - pi/2. g(0) = -pi/2. g'(x)=tan^2(x)>=0. So g is increasing.
    # g(pi/2 - eps) = cot(eps) - (pi/2-eps) - pi/2 = 1/eps - pi + eps. This is positive for small eps.
    # So there is a unique root in (0, pi/2). My first solver call was right after all.

    sol = root(func, 1.2) # Initial guess x=1.2
    x0 = sol.x[0]

    # The constant C is cot^2(x0/2)
    C = (1 / np.tan(x0 / 2))**2
    
    print("To find the constant C, we solve the transcendental equation for a value x0:")
    print("tan(x0) = x0 + pi/2")
    print(f"The numerical solution is x0 = {x0:.8f}")
    print("\nThe constant C is then calculated using the formula C = cot(x0/2)^2:")
    print(f"C = cot({x0:.8f} / 2)^2")
    print(f"C = {C:.8f}")

if __name__ == '__main__':
    get_constant_c()