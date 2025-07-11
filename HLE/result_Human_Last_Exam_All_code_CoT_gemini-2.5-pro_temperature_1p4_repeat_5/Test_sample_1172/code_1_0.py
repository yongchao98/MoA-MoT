import sympy as sp

def solve_mutual_inductance_change():
    """
    This function calculates the change in mutual inductance per unit length
    between two circuits when a cylindrical magnetic shield is added.

    The final expression for the change in mutual inductance per unit length (Delta_M') is derived and printed.
    """

    # Define the symbols for our physical quantities.
    # mu_0: Permeability of free space
    # h: Separation of wires within a circuit
    # d: Distance between the centers of the two circuits
    # R1: Inner radius of the cylindrical shell
    mu_0, h, d, R1 = sp.symbols('mu_0 h d R1', real=True, positive=True)
    pi = sp.pi

    # --- Step 1: Model the circuits as magnetic dipoles ---
    # Circuit 1 is at x = -d/2, Circuit 2 is at x = d/2.
    # Both circuits have currents +I and -I, with wire separation h.
    # Let's assume the wire with +I is to the left of the wire with -I in both circuits.
    # The magnetic dipole moment vector points from the -I wire to the +I wire.
    # So, for Circuit 1 (+I at -d/2-h/2, -I at -d/2+h/2), the vector is from right to left.
    # m1 = -I*h*x_hat
    # Similarly for Circuit 2, m2 = -I*h*x_hat

    # --- Step 2: Model the effect of the shield using the method of images ---
    # The shield with the given properties acts as a perfect magnetic shield,
    # which can be modeled by placing an image dipole (m_image) to satisfy the
    # boundary condition that the magnetic potential is zero at r=R1.

    # Position of image dipole: The image of a source at a from the center of a
    # cylinder of radius R1 is at b = R1^2/a.
    # Circuit 1 is at a = d/2 from the center.
    # The image is at b = R1**2 / (d/2) = 2*R1**2 / d, on the same line from the origin.
    # So, r_image = (-2*R1**2/d, 0).

    # Moment of the image dipole: For two wires, the image dipole moment is found to be:
    # m_image = -I * h * (4*R1**2 / d**2) * x_hat.
    # Let's represent the moments divided by current I and the vector part.
    m_image_moment_mag = h * (4 * R1**2 / d**2)
    m2_moment_mag = h
    
    # --- Step 3: Calculate the mutual inductance between the image and circuit 2 ---
    # The change in mutual inductance is the inductance between m_image and m2.
    # The formula for mutual inductance per unit length between two dipoles m_a, m_b is:
    # M' = (mu_0 / (2*pi*r^2)) * ( (m_a.m_b)/I^2 - 2*(m_a.r_hat)*(m_b.r_hat)/I^2 )
    
    # The vector from m_image to m2 is r_vec:
    # r_vec = (d/2 - (-2*R1**2/d), 0) = (d/2 + 2*R1**2/d, 0)
    # The magnitude squared of this vector is r_sq:
    r_sq = (d/2 + 2*R1**2/d)**2

    # The unit vector r_hat is in the +x direction.
    # The dipole moments are in the -x direction.

    # Calculate the dot products needed (normalized by I^2):
    # (m_image . m2)/I^2 = (-m_image_mag * x_hat) . (-m2_mag * x_hat)
    m_image_dot_m2_norm = m_image_moment_mag * m2_moment_mag
    
    # (m_image . r_hat)/I = (-m_image_mag * x_hat) . (x_hat)
    m_image_dot_r_hat_norm = -m_image_moment_mag
    
    # (m2 . r_hat)/I = (-m2_mag * x_hat) . (x_hat)
    m2_dot_r_hat_norm = -m2_moment_mag
    
    # Combine them for the interaction term in the formula
    interaction_term = m_image_dot_m2_norm - 2 * m_image_dot_r_hat_norm * m2_dot_r_hat_norm
    # interaction_term = (m_image_mag * m2_mag) - 2 * (-m_image_mag) * (-m2_mag)
    #                  = m_image_mag * m2_mag - 2 * m_image_mag * m2_mag
    #                  = -m_image_mag * m2_mag
    interaction_term_val = - (h * (4*R1**2 / d**2)) * h
    
    # Assemble the final expression for the change in mutual inductance per unit length
    Delta_M_prime = (mu_0 / (2 * pi * r_sq)) * interaction_term_val

    # --- Step 4: Simplify the expression and print ---
    final_expression = sp.simplify(Delta_M_prime)

    # Print the final result in a human-readable format
    # The problem asks to output each number, so we construct the string carefully.
    numerator = sp.latex( -final_expression.args[0] * final_expression.args[1])
    denominator = sp.latex(final_expression.args[2]**-1)
    
    print("The change in mutual inductance per unit length, Delta_M', is:")
    
    # A simplified and readable printout
    print("\nDelta_M' = - (8 * mu_0 * h**2 * R1**2) / (pi * (d**2 + 4*R1**2)**2)")
    
    # For verification, print the symbolic object as well
    # print("\nSympy object:")
    # sp.pprint(final_expression)

if __name__ == '__main__':
    solve_mutual_inductance_change()
    # The final expression represents the change in mutual inductance per unit length.
    final_answer_expression = "-8 * mu_0 * h**2 * R1**2 / (pi * (d**2 + 4 * R1**2)**2)"
    # The format specification is <<<answer content>>>
    # Since the answer is an expression, we will put it in that format.
    print(f"\n<<<{final_answer_expression}>>>")
