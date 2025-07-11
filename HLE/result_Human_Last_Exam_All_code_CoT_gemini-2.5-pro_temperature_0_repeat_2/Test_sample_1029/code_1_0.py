def solve_poynting_vector():
    """
    This function derives and prints the Poynting vector for a moving,
    charged cylindrical rod in an external electric field.
    """
    # Use unicode characters for better readability
    rho = "\u03C1"  # ρ
    mu0 = "\u03BC\u2080"  # μ₀
    eps0 = "\u03B5\u2080"  # ε₀
    phi_hat = "\u03C6\u0302" # φ̂
    r_hat = "r\u0302"   # r̂
    k_hat = "k\u0302"   # k̂
    vec_S = "\u20D7S"   # Vector S
    vec_E = "\u20D7E"   # Vector E
    vec_B = "\u20D7B"   # Vector B
    vec_J = "\u20D7J"   # Vector J
    vec_v = "\u20D7v"   # Vector v

    print("### Computation of the Poynting Vector ###\n")
    print(f"The Poynting vector {vec_S} is defined as: {vec_S} = (1/{mu0}) * ({vec_E} x {vec_B})\n")
    print("We will solve this step-by-step in a cylindrical coordinate system (r, \u03C6, z).\n")

    # --- Step 1: Total Electric Field ---
    print("--- Step 1: Finding the Total Electric Field " + vec_E + " ---")
    print(f"The total electric field is the sum of the external field and the rod's own field: {vec_E} = {vec_E}_ext + {vec_E}_rod.")
    print(f"The external field is given as uniform and along the axis: {vec_E}_ext = E {k_hat}.")
    print(f"The rod's field {vec_E}_rod is found using Gauss's Law. It points radially.")
    print(f"  - Inside the rod (r <= R): {vec_E}_rod = ({rho}*r / (2*{eps0})) {r_hat}")
    print(f"  - Outside the rod (r > R): {vec_E}_rod = ({rho}*R\u00B2 / (2*{eps0}*r)) {r_hat}")
    print("Therefore, the total electric field is:")
    print(f"  - Inside (r <= R):   {vec_E} = E {k_hat} + ({rho}*r / (2*{eps0})) {r_hat}")
    print(f"  - Outside (r > R):  {vec_E} = E {k_hat} + ({rho}*R\u00B2 / (2*{eps0}*r)) {r_hat}\n")

    # --- Step 2: Magnetic Field ---
    print("--- Step 2: Finding the Magnetic Field " + vec_B + " ---")
    print(f"The rod moves with velocity {vec_v} = v {k_hat}, creating a volume current density {vec_J} = {rho}*{vec_v} = {rho}*v {k_hat}.")
    print(f"This current generates a magnetic field {vec_B}, found using Ampere's Law. The field lines are circles around the axis ({phi_hat} direction).")
    print(f"  - Inside the rod (r <= R): {vec_B} = ({mu0}*{rho}*v*r / 2) {phi_hat}")
    print(f"  - Outside the rod (r > R): {vec_B} = ({mu0}*{rho}*v*R\u00B2 / (2*r)) {phi_hat}\n")

    # --- Step 3: Poynting Vector ---
    print("--- Step 3: Computing the Poynting Vector " + vec_S + " ---")
    print(f"Now we compute the cross product {vec_E} x {vec_B} and divide by {mu0}.")
    print(f"We use the cross products: {k_hat} x {phi_hat} = -{r_hat} and {r_hat} x {phi_hat} = {k_hat}.\n")

    print(">>> Final Result <<<\n")
    # For r <= R
    print("1. Inside the rod (r <= R):")
    print(f"{vec_S}_in = (1/{mu0}) * [ (E {k_hat} + ({rho}*r / (2*{eps0})) {r_hat}) x (({mu0}*{rho}*v*r / 2) {phi_hat}) ]")
    print(f"{vec_S}_in = (1/{mu0}) * [ E*({mu0}*{rho}*v*r/2)({k_hat} x {phi_hat}) + ({rho}*r/(2*{eps0}))*({mu0}*{rho}*v*r/2)({r_hat} x {phi_hat}) ]")
    print(f"{vec_S}_in = (E*{rho}*v*r/2)(- {r_hat}) + ({rho}\u00B2*v*r\u00B2/(4*{eps0}))({k_hat})")
    print("\nFinal Equation (Inside):")
    print(f"S_radial = - (E * {rho} * v * r / 2)")
    print(f"S_axial  =   ({rho}\u00B2 * v * r\u00B2 / (4 * {eps0}))")
    print(f"{vec_S}_in = - (E*{rho}vr/2) {r_hat} + ({rho}\u00B2vr\u00B2/(4{eps0})) {k_hat}\n")

    # For r > R
    print("2. Outside the rod (r > R):")
    print(f"{vec_S}_out = (1/{mu0}) * [ (E {k_hat} + ({rho}*R\u00B2 / (2*{eps0}*r)) {r_hat}) x (({mu0}*{rho}*v*R\u00B2 / (2*r)) {phi_hat}) ]")
    print(f"{vec_S}_out = (1/{mu0}) * [ E*({mu0}*{rho}*v*R\u00B2/(2r))({k_hat} x {phi_hat}) + ({rho}*R\u00B2/(2*{eps0}r))*({mu0}*{rho}*v*R\u00B2/(2r))({r_hat} x {phi_hat}) ]")
    print(f"{vec_S}_out = (E*{rho}*v*R\u00B2/(2r))(- {r_hat}) + ({rho}\u00B2*v*R\u2074/(4*{eps0}r\u00B2))({k_hat})")
    print("\nFinal Equation (Outside):")
    print(f"S_radial = - (E * {rho} * v * R\u00B2 / (2 * r))")
    print(f"S_axial  =   ({rho}\u00B2 * v * R\u2074 / (4 * {eps0} * r\u00B2))")
    print(f"{vec_S}_out = - (E*{rho}vR\u00B2/(2r)) {r_hat} + ({rho}\u00B2vR\u2074/(4{eps0}r\u00B2)) {k_hat}\n")

if __name__ == '__main__':
    solve_poynting_vector()