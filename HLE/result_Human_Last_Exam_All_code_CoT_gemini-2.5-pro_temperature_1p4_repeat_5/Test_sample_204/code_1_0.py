import numpy as np

def solve_hopf_charge():
    """
    Calculates the Hopf charge of the given vector field by finding the linking
    number of the preimages of two regular points.
    """
    
    # The vector field (nx,ny,nz) is defined by:
    # f = atan2(y,x)
    # r2 = sqrt((x*x+y*y-0.5)**2 + z*z)
    # G = PI * exp(-10*r2)
    # n = (sin(G)*cos(f), sin(G)*sin(f), cos(G))
    # where PI is the mathematical constant pi.
    
    print("Step 1: The Hopf charge H is the linking number of the preimages of two regular values.")
    
    # Step 2: Find the preimage C1 for the point p1 = (0,0,-1).
    # n = (0,0,-1) implies cos(G) = -1, which means G = pi.
    # G = pi * exp(-10*r2) = pi  =>  exp(-10*r2) = 1  =>  r2 = 0.
    # r2 = sqrt((x*x+y*y-0.5)**2 + z*z) = 0
    # This implies z=0 and (x*x+y*y-0.5) = 0.
    # So, C1 is the circle: x^2 + y^2 = 0.5, z = 0.
    print("Step 2: The preimage C1 for n=(0,0,-1) is the circle x^2+y^2=0.5 in the z=0 plane.")
    
    # Step 3: Find the preimage C2 for the point p2 = (1,0,0).
    # n = (1,0,0) implies cos(G) = 0 and (cos(f)=1, sin(f)=0).
    # f = 0 implies y = 0 and x > 0.
    # cos(G) = 0 implies G = pi/2.
    # G = pi * exp(-10*r2) = pi/2  =>  exp(-10*r2) = 0.5  =>  r2 = ln(2)/10.
    # With y=0, r2^2 = (x^2 - 0.5)^2 + z^2.
    # So, C2 is the curve: (x^2 - 0.5)^2 + z^2 = (ln(2)/10)^2, for y=0 and x>0.
    R_c2_sq = (np.log(2)/10)**2
    print(f"Step 3: The preimage C2 for n=(1,0,0) is the loop (x^2-0.5)^2 + z^2 = {R_c2_sq:.6f}, in the y=0 plane.")

    # Step 4: Determine the linking number.
    # Curve C1 is a circle of radius sqrt(0.5) in the xy-plane. It pierces the xz-plane (which contains C2) at the point (sqrt(0.5), 0, 0).
    # Curve C2 is a loop in the xz-plane. The equation (x^2 - 0.5)^2 + z^2 = R^2 shows that this loop is centered around x^2 = 0.5 and z=0.
    # This means C2 is a loop that encircles the point where C1 pierces its plane.
    # Two such curves are linked exactly once.
    hopf_charge = 1
    print("Step 4: The curve C2 forms a closed loop around the point where C1 passes through its plane. Their linking number is 1.")

    print("\nFinal Result Equation:")
    print(f"Hopf Charge = {hopf_charge}")

solve_hopf_charge()