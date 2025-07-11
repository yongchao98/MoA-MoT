import numpy as np

def calculate_shape():
    """
    This function demonstrates the shape of the set S for a 2D example.
    It generates points in S and verifies they lie on an ellipse.
    """
    # Define the angle theta for the second vector y2
    theta = np.pi / 4

    # The set S is parameterized by an angle phi for the unit vector s
    # S = { (cos^2(phi), cos^2(phi - theta)) | phi in [0, 2pi) }
    
    # We derived the implicit equation for this set:
    # (x + y - sin^2(theta))^2 = 4*x*y*cos^2(theta)
    # For theta = pi/4, this simplifies to the circle:
    # (x - 0.5)^2 + (y - 0.5)^2 = 0.25

    print("For the 2D example with y1=(1,0) and y2=(cos(pi/4), sin(pi/4)):")
    print("The set S is described by the parametric equations:")
    print("x = cos(phi)^2")
    print("y = cos(phi - pi/4)^2")
    print("These points lie on a circle (a special type of ellipse) defined by the equation:")
    print("(x - 0.5)^2 + (y - 0.5)^2 = 0.25")
    
    # Expanding the equation: x^2 - x + y^2 - y + 0.25 = 0
    
    # Per instructions, outputting each number in the final equation.
    # Using the form (x - a)^2 + (y - b)^2 = r^2
    a = 0.5
    b = 0.5
    r_squared = 0.25
    print("\nThe numbers in the equation (x - a)^2 + (y - b)^2 = r^2 are:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"r^2 = {r_squared}")

    # Verification for a few points
    print("\nVerifying for a few points:")
    for phi_deg in [0, 45, 90, 135]:
        phi = np.deg2rad(phi_deg)
        x = np.cos(phi)**2
        y = np.cos(phi - theta)**2
        
        # Check if the point satisfies the circle equation
        lhs = (x - 0.5)**2 + (y - 0.5)**2
        rhs = 0.25
        
        print(f"For phi = {phi_deg} degrees:")
        print(f"  x = {x:.4f}, y = {y:.4f}")
        print(f"  (x-0.5)^2 + (y-0.5)^2 = {lhs:.4f} (should be close to {rhs})")

calculate_shape()
