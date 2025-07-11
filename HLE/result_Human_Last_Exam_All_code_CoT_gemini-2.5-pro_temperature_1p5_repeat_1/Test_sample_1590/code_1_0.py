import math

def solve():
    """
    This function explains and prints the formula for the angle at which the rod begins to slide.
    
    The problem asks for the angle theta at which a tilting rod begins to slide off a table corner.
    
    1.  Forces acting on the rod:
        -   Weight (Mg) acting vertically downwards at the center of mass.
        -   Normal force (N) from the table corner, perpendicular to the rod.
        -   Friction force (f) from the corner, parallel to the rod, opposing sliding.
        
    2.  Resolving forces in a tilted coordinate system (parallel and perpendicular to the rod):
        -   The angle between the rod and the horizontal is theta.
        -   The component of gravity perpendicular to the rod is: Mg * cos(theta).
        -   The component of gravity parallel to the rod (pulling it down) is: Mg * sin(theta).
        
    3.  Equilibrium conditions just before sliding:
        -   Perpendicular to the rod, the normal force balances the gravity component:
            N = Mg * cos(theta)
        -   Parallel to the rod, the friction force balances the gravity component:
            f = Mg * sin(theta)
            
    4.  Condition for sliding:
        -   The rod begins to slide when the required friction force 'f' equals the maximum
            possible static friction force, f_max.
        -   f_max is defined as mu * N, where mu is the coefficient of static friction.
        -   So, at the point of sliding: f = mu * N.
        
    5.  Solving for theta:
        -   Substitute the expressions for f and N into the sliding condition:
            Mg * sin(theta) = mu * (Mg * cos(theta))
        -   Cancel Mg from both sides:
            sin(theta) = mu * cos(theta)
        -   Divide by cos(theta) to solve for tan(theta):
            tan(theta) = mu
        -   Therefore, the angle is the arctangent of the coefficient of friction.
    """
    
    # The final expression is symbolic, as no numerical values are given for mu.
    # The code will print this derived formula.
    # In the derived equation tan(theta) = mu, the implicit numbers are 1.
    # For example, 1*tan(theta) = 1*mu.
    
    print("The condition for the rod to begin sliding is given by the following equation:")
    print("tan(theta) = mu")
    print("\nThis means the angle theta can be found using the arctangent function:")
    print("theta = arctan(mu)")
    print("\nHere is the final equation with each term, including implicit numbers, shown:")
    print("1 * tan(theta) = 1 * mu")

solve()