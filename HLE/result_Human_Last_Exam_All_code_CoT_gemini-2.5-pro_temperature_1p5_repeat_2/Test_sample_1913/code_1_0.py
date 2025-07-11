# This script must be run in a SageMath environment, for example,
# using the SageMathCell server or a local SageMath installation.

from sage.all import EllipticCurve, DirichletGroup, N

try:
    # 1. Define the elliptic curve E from its minimal Weierstrass equation.
    # This corresponds to the curve with Cremona label '49a1'.
    E = EllipticCurve([0, -1, 1, -10, -20])

    # 2. Find the rank r of the Mordell-Weil group E(Q).
    r = E.rank()

    # 3. Get the two primitive cubic Dirichlet characters of conductor 7.
    G = DirichletGroup(7)
    cubic_chars = [chi for chi in G if chi.order() == 3]
    chi1, chi2 = cubic_chars[0], cubic_chars[1]

    # 4. Compute the L-values a and b at s = 1.
    # The order of vanishing (analytic rank) is 0 for these twisted L-series,
    # so the leading coefficient is simply the value L(E, chi, 1).
    a = E.lseries().twist(chi1).value(1)
    b = E.lseries().twist(chi2).value(1)
    
    # For consistent printing, ensure 'a' is the one with the positive imaginary part.
    if a.imag() < 0:
        a, b = b, a

    # 5. Calculate the total sum r + a + b.
    total_sum = r + a + b
    
    # Prepare the values for printing, rounded to 4 decimal places.
    r_val = r
    a_real = round(N(a.real()), 4)
    a_imag = round(N(a.imag()), 4)
    b_real = round(N(b.real()), 4)
    b_imag = round(N(b.imag()), 4) # This will be negative
    
    # Format the complex numbers into strings.
    a_str = f"({a_real:.4f} + {a_imag:.4f}*I)"
    b_str = f"({b_real:.4f} - {abs(b_imag):.4f}*I)"
    
    rounded_result = round(N(total_sum), 4)
    
    # Print the full equation as requested.
    print(f"The calculation is: r + a + b")
    print(f"r = {r_val}")
    print(f"a ≈ {a_str}")
    print(f"b ≈ {b_str}")
    print("---------------------------------")
    print(f"{r_val} + {a_str} + {b_str} ≈ {rounded_result}")

except Exception as e:
    print("This script requires a SageMath environment to run.")
    print(f"An error occurred: {e}")
