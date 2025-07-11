import numpy as np

def solve():
    """
    This function determines the function h(x) and prints its coefficients.
    
    The condition for a(t) -> 0 is that the initial point (a(0), b(0)) lies inside the separatrix
    that bounds the basin of attraction of the stable node (0, 1/2). We are looking for the equation
    of this separatrix, a^2 = h(b).

    A curve a^2 = h(b) is a trajectory if d/dt(a^2 - h(b)) = 0.
    d/dt(a^2 - h(b)) = 2*a*a'(t) - h'(b)*b'(t)
                     = 2*a*(-1/2*a^2 - 2*b^2 - b + 1) - h'(b)*(-a*b)
                     = -a^3 - 4*a*b^2 - 2*a*b + 2*a + a*b*h'(b)
                     = a * (-a^2 - 4*b^2 - 2*b + 2 + b*h'(b))

    Setting the term in the parenthesis to zero, and substituting a^2 = h(b):
    -h(b) - 4*b^2 - 2*b + 2 + b*h'(b) = 0
    This gives a linear first-order ODE for h(b):
    b*h'(b) - h(b) = 4*b^2 + 2*b - 2

    This ODE has a general solution of the form h(b) = h_p(b) + h_h(b), where h_p is a particular solution
    and h_h is the solution to the homogeneous equation b*h' - h = 0.
    The homogeneous solution is h_h(b) = C*b for some constant C.
    A particular solution can be found to be h_p(b) = 4*b^2 + 2*b*ln(b) + 2.
    So, the general solution is h(b) = 4*b^2 + 2*b*ln(b) + 2 + C*b.

    The separatrix must connect the saddle points (±sqrt(2), 0). This means that as b -> 0,
    a^2 must approach 2. So, we need lim_{b->0} h(b) = 2.
    Let's check our general solution:
    lim_{b->0} (4*b^2 + 2*b*ln(b) + 2 + C*b) = 0 + 0 + 2 + 0 = 2.
    This condition is satisfied for any constant C. Since no other information is provided to
    determine C, we choose the simplest and most natural option, C=0.

    Thus, the function is h(x) = 4x^2 + 2x*ln(x) + 2.
    """
    
    # Coefficients for the function h(x) = c1*x^2 + c2*x*ln(x) + c3*x + c4
    c1 = 4.0
    c2 = 2.0
    c3 = 0.0
    c4 = 2.0

    print("The function h(x) is determined to be of the form:")
    print("h(x) = c1*x^2 + c2*x*ln(x) + c3*x + c4")
    print("\nThe coefficients are:")
    print(f"c1 = {c1}")
    print(f"c2 = {c2}")
    print(f"c3 = {c3}")
    print(f"c4 = {c4}")
    
    # Final answer format requires showing the final equation with numbers
    print("\nSo the final equation for h(x) is:")
    print(f"h(x) = {c1}*x^2 + {c2}*x*ln(x) + {c4}")

solve()