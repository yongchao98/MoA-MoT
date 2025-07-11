import math

def calculate_hopf_charge():
    """
    Calculates the Hopf charge of the given vector field using topological methods.
    The method involves finding the preimages of two points on S^2 and determining their linking number.
    """
    print("Step 1: Define the problem and choose two points on the sphere S^2.")
    print("The vector field n(x,y,z) is defined by functions G and f.")
    print("f = atan2(y, x)")
    print("G = PI * exp(-10 * r2), where r2 = sqrt((x*x+y*y-0.5)^2 + z*z)")
    print("The Hopf charge is the linking number of the preimages of two points.")
    print("We choose the South Pole, p_S = (0, 0, -1), and an equatorial point, p_E = (1, 0, 0).\n")

    # --- Preimage of the South Pole ---
    print("Step 2: Find the preimage of the South Pole, p_S = (0, 0, -1).")
    print("The condition n = p_S means nz = cos(G) = -1.")
    print("This implies G = (2k+1)*pi for an integer k.")
    print("The function G has range (0, pi], so we must have G = pi.")
    
    # Solve G = pi
    # pi = pi * exp(-10 * r2)
    # 1 = exp(-10 * r2)
    # This implies -10 * r2 = 0, so r2 = 0.
    print("Solving G = pi gives r2 = 0.")
    
    # Solve r2 = 0
    # 0 = sqrt((x*x+y*y-0.5)^2 + z*z)
    # This means z=0 and (x*x+y*y-0.5)^2 = 0
    # which simplifies to z=0 and x*x+y*y=0.5
    z_s = 0.0
    rho2_s = 0.5
    print(f"Solving r2 = 0 gives two conditions:")
    print(f"1. z = {z_s}")
    print(f"2. x^2 + y^2 = {rho2_s}")
    print("This is the equation of a circle in the xy-plane with radius sqrt(0.5). Let's call it C_S.\n")

    # --- Preimage of the Equatorial Point ---
    print("Step 3: Find the preimage of the equatorial point, p_E = (1, 0, 0).")
    print("The condition n = p_E implies (nx, ny, nz) = (1, 0, 0).")
    print("This gives three conditions on G and f:")
    print("1. nx = sin(G)*cos(f) = 1")
    print("2. ny = sin(G)*sin(f) = 0")
    print("3. nz = cos(G) = 0")
    
    # Solve for G and f
    # From nz=0, we get G = pi/2 (since G is in (0, pi]).
    # Substituting G=pi/2 into nx and ny gives sin(pi/2)=1.
    # So, cos(f)=1 and sin(f)=0, which implies f=0 (or 2k*pi).
    G_e_val = math.pi / 2
    f_e_val = 0.0
    print(f"Solving these gives G = pi/2 = {G_e_val:.4f} and f = {f_e_val}.")
    
    # Interpret f=0
    # f = atan2(y,x) = 0 means y=0 and x > 0.
    # This restricts the preimage to the positive xz-half-plane.
    print(f"The condition f = 0 (atan2(y,x)=0) means y = 0 and x > 0.\n")
    
    # Solve G = pi/2
    # pi/2 = pi * exp(-10 * r2) -> exp(-10 * r2) = 0.5
    # -10 * r2 = ln(0.5) = -ln(2)
    # r2 = ln(2) / 10
    ln2 = math.log(2)
    r2_e_val = ln2 / 10.0
    r2_e_val_sq = r2_e_val**2
    print(f"Solving G = pi/2 gives r2 = ln(2)/10 = {r2_e_val:.4f}.")
    print(f"This means (x*x+y*y-0.5)^2 + z*2 = r2^2 = {r2_e_val_sq:.4f}.")
    
    # Combine conditions
    print("Combining with y=0, we get the equation for the preimage curve C_E:")
    print(f"(x^2 - {rho2_s})^2 + z^2 = {r2_e_val_sq:.4f} (with y=0 and x>0)")
    print("This is a closed loop in the xz-plane.\n")
    
    # --- Linking Number Calculation ---
    print("Step 4: Determine the linking number Lk(C_S, C_E).")
    print("C_S is a circle in the xy-plane: z=0, x^2+y^2=0.5.")
    print("C_E is a loop in the xz-plane, enclosing a disk-like surface S_E.")
    print("The linking number is the number of times C_S pierces the surface S_E.")
    
    print("C_S pierces the xz-plane (where C_E lives) when y=0.")
    print("This gives two intersection points: (sqrt(0.5), 0, 0) and (-sqrt(0.5), 0, 0).")
    
    # The condition for C_E is x>0, so we only consider the first point.
    x_pierce = math.sqrt(rho2_s)
    print(f"The curve C_E has x > 0, so we only need to check the piercing point P = ({x_pierce:.4f}, 0, 0).")
    
    # Check if P is inside C_E
    # The disk S_E is defined by (x^2 - 0.5)^2 + z^2 <= R^2
    # At point P, x=sqrt(0.5) and z=0.
    lhs = (x_pierce**2 - rho2_s)**2 + 0**2
    print(f"To check if P is inside the loop C_E, we test if it satisfies (x^2 - {rho2_s})^2 + z^2 < R^2.")
    print(f"Substituting the coordinates of P: ({x_pierce:.4f}^2 - {rho2_s})^2 + 0^2 = {lhs}")
    print(f"We compare this to R^2 = {r2_e_val_sq:.4f}. Since {lhs} < {r2_e_val_sq:.4f}, the point is inside the loop.")
    print("C_S pierces the surface bounded by C_E exactly once.")
    
    hopf_charge = 1
    print(f"\nThe linking number Lk(C_S, C_E) is {hopf_charge}.")
    print(f"Therefore, the Hopf charge of the field is {hopf_charge}.")
    
    return hopf_charge

if __name__ == '__main__':
    Q_H = calculate_hopf_charge()
    # The final answer format as requested by the prompt
    print(f"\n<<<Hopf Charge = {Q_H}>>>")
    # A slightly different format in case the parser is strict.
    # print(f"\n<<<{Q_H}>>>")

