import textwrap

def solve_toroidal_expansion():
    """
    Identifies and explains the correct spectral series for poloidal dependence in toroidal systems.
    """
    options = {
        'A': 'Gegenbauer Polynomials',
        'B': 'Spherical harmonics expansion',
        'C': 'B-splines',
        'D': 'Fourier series',
        'E': 'Chebyshev polynomials',
        'F': 'Spherical Harmonics', # Duplicate concept of B
        'G': 'Hermite Polynomials',
        'H': 'Jacobi Polynomials',
        'I': 'Legendre Polynomials',
        'J': 'Fourier-Legendre Series'
    }

    correct_option_key = 'D'
    correct_option_value = options[correct_option_key]

    explanation = """
    In toroidal systems, such as tokamaks in fusion research, physical quantities are described using a coordinate system where one angle, the poloidal angle (θ), represents the position 'the short way around' the torus.

    This poloidal angle is periodic, meaning the system returns to the same state after a rotation of 2π radians. The most natural and widely-used mathematical technique to represent any periodic function is a Fourier series expansion.

    A function f(θ) dependent on the poloidal angle is thus expanded as a sum of sines and cosines (or complex exponentials).
    A simplified form of the expansion equation is:
    f(θ) = Σ [ A_m * cos(m*θ) + B_m * sin(m*θ) ]
    where 'm' is the poloidal mode number.

    Therefore, the correct spectral expansion technique is the Fourier series.
    """

    print(textwrap.dedent(explanation).strip())
    print("\n----------------------------------")
    print(f"The correct option is: {correct_option_key}. {correct_option_value}")

if __name__ == "__main__":
    solve_toroidal_expansion()