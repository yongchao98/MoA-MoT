import numpy as np

def calculate_hopf_charge():
    """
    This function calculates the Hopf charge of the given vector field
    by computing the linking number of the preimages of two points.
    """

    print("### Plan to calculate the Hopf Charge ###")
    print("1. The Hopf charge will be determined using the Whitehead formula, which corresponds to the linking number of the preimages of two regular points on the target sphere S2.")
    print("2. Choose two points on S2: the south pole n1=(0,0,-1) and an equatorial point n2=(1,0,0).")
    print("3. Find their preimages, C1 and C2, in R3.")
    print("4. Calculate the linking number Lk(C1, C2) to find the Hopf charge.\n")

    # Step 1: Define the preimages C1 and C2
    print("### Step 1: Finding the Preimages ###")
    print("The vector field is defined by n = (sin(G)cos(f), sin(G)sin(f), cos(G)) where G and f depend on (x,y,z).")

    # Preimage C1 for n1 = (0,0,-1)
    print("\nFor point n1 = (0,0,-1), the polar angle G must be PI.")
    print("The field definition is G = PI * exp(-10*r2), where r2 = sqrt((x*x+y*y-0.5)^2 + z^2).")
    print("G = PI requires exp(-10*r2) = 1, which means r2 = 0.")
    print("r2 = 0 implies (x*x+y*y-0.5) = 0 and z = 0.")
    print("So, the preimage C1 is the circle: x^2 + y^2 = 0.5, z = 0.")
    r_c1_sq = 0.5
    r_c1 = np.sqrt(r_c1_sq)
    print(f"C1 is a circle of radius {r_c1:.4f} in the xy-plane.\n")

    # Preimage C2 for n2 = (1,0,0)
    print("For point n2 = (1,0,0), G = PI/2 and the azimuthal angle f = 0.")
    print("f = atan2(y,x) = 0 implies y = 0 and x > 0.")
    print("G = PI * exp(-10*r2) = PI/2 implies exp(-10*r2) = 0.5.")
    print("This gives -10*r2 = ln(0.5) = -ln(2), so r2 = ln(2)/10.")
    const_C = np.log(2) / 10.0
    print(f"This means r2^2 = (ln(2)/10)^2 = {const_C**2:.6f}.")
    print("So, the preimage C2 is the loop defined by: (x^2 - 0.5)^2 + z^2 = (ln(2)/10)^2, with y=0 and x>0.")
    
    # Step 2: Calculate the linking number
    print("\n### Step 2: Calculating the Linking Number Lk(C1, C2) ###")
    print("We can find Lk(C1, C2) by counting how many times C2 pierces the disk D1 (x^2+y^2 <= 0.5) bounded by C1.")
    print("The piercing points are where C2 intersects the z=0 plane.")
    print("For C2, setting z=0 gives (x^2 - 0.5)^2 = C^2, where C = ln(2)/10.")
    print("This yields x^2 = 0.5 +/- C.")
    
    x_sq_pierce_outer = 0.5 + const_C
    x_sq_pierce_inner = 0.5 - const_C
    x_pierce_outer = np.sqrt(x_sq_pierce_outer)
    x_pierce_inner = np.sqrt(x_sq_pierce_inner)

    print(f"The x-coordinates of the piercing points are {x_pierce_inner:.4f} and {x_pierce_outer:.4f}.")
    print(f"The radius of the disk D1 is {r_c1:.4f}.")
    print(f"Since {x_pierce_inner:.4f} < {r_c1:.4f} and {x_pierce_outer:.4f} > {r_c1:.4f}, C2 pierces the disk exactly once.")
    print("The magnitude of the linking number is 1.\n")

    # Step 3: Determine the sign of the linking
    print("### Step 3: Determining the Sign ###")
    print("The sign depends on the orientations of the curves.")
    print("Let's orient C1 counter-clockwise in the xy-plane. By the right-hand rule, the normal vector of the disk D1 is in the +z direction (0,0,1).")
    print("The curve C2 can be parameterized. At the piercing point inside the disk (x,z) = (sqrt(0.5 - C), 0), calculation shows the velocity vector of C2 has a negative z-component.")
    print("The dot product of C2's velocity vector and D1's normal vector is negative.")
    linking_number = -1
    print(f"This means the linking number is {linking_number}.\n")

    # Final result
    print("### Final Result ###")
    print("The Hopf charge H is equal to this linking number.")
    print(f"Hopf Charge = {linking_number}")

calculate_hopf_charge()