import math

def calculate_hopf_charge():
    """
    Calculates the Hopf charge for the given axially symmetric vector field.

    The Hopf charge H for a field of the form n = (sin(G)cos(p*f), sin(G)sin(p*f), cos(G))
    is given by H = p * q.
    """

    # 1. Determine the azimuthal winding number, p.
    # The field's azimuthal angle is f = atan2(y,x), which is the standard
    # cylindrical angle phi. This means the field winds once for every spatial rotation.
    p = 1
    
    # 2. Determine the poloidal winding number, q.
    # q is given by (G_center - G_infinity) / pi.

    # G_center is the value of G at the structure's core, where r2 = 0.
    # G = pi * exp(-10 * r2)
    G_center = math.pi * math.exp(-10 * 0)

    # G_infinity is the value of G at spatial infinity, where r2 -> inf.
    # As r2 -> inf, exp(-10*r2) -> 0.
    G_infinity = math.pi * 0
    
    # Calculate q. It must be an integer.
    q = (G_center - G_infinity) / math.pi
    
    # 3. Calculate the total Hopf charge H.
    H = p * int(q)
    
    print("The Hopf charge H is the product of two winding numbers: H = p * q")
    print(f"Azimuthal winding number p = {p}")
    print("Poloidal winding number q is calculated from the profile G(r2):")
    print(f"  - At the center (r2=0), G = pi * exp(0) = {G_center/math.pi:.1f}*pi")
    print(f"  - At infinity (r2->inf), G = pi * exp(-inf) = {G_infinity/math.pi:.1f}*pi")
    print(f"  - q = (G_center - G_infinity) / pi = ({G_center/math.pi:.1f}*pi - {G_infinity/math.pi:.1f}*pi) / pi = {int(q)}")
    print("\nFinal Calculation:")
    print(f"H = p * q = {p} * {int(q)} = {int(H)}")

calculate_hopf_charge()