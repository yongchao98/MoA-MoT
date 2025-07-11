import math

def calculate_time_dilation():
    """
    This function calculates the gravitational time dilation factor f
    and the C program's variable memory usage z, following the constraints
    of the Wuxing architecture.
    """

    # 1. Approximate physical constants using the frac format (n/d * 10**e).
    #    n and d must fit in 2 decimal digits.
    # G = 6.674e-11 m^3 kg^-1 s^-2. Approx: 60/9 * 10**-11
    g_val = (60 / 9) * (10**-11)

    # M_pulsar = 2 * M_sun = 3.978e30 kg. Approx: 79/20 * 10**30
    m_val = (79 / 20) * (10**30)

    # r = Pandora_radius + distance = 20km + 60km = 80km = 8e4 m.
    r_val = 8 * (10**4)

    # c = 2.998e8 m/s. Approx: 3 * 10**8
    c_val = 3 * (10**8)

    # 2. Calculate the term x = (2*G*M)/(r*c^2).
    #    The C program would store this intermediate result in a register
    #    as part of the calculation for 'A'.
    x = (2 * g_val * m_val) / (r_val * c_val**2)

    # 3. Calculate A = 1 - x. This would be stored in the first 'frac' variable.
    A = 1 - x

    # 4. Use Newton's method to find sqrt(A), since sqrt() is unavailable.
    #    y_{n+1} = (y_n + A / y_n) / 2.
    #    The second 'frac' variable would hold the result 'f'.
    #    Three iterations provide sufficient precision.
    f = 1.0  # Initial guess for the square root.
    f = (f + A / f) / 2.0
    f = (f + A / f) / 2.0
    f = (f + A / f) / 2.0

    # Round final result to 0.0001 precision.
    f_rounded = round(f, 4)

    # 5. Calculate the memory usage 'z' for variables in the C program.
    #    The most efficient C program uses two 'frac' variables.
    #    A 'frac' contains three 'char' types, each 2 decimal digits (2D).
    #    Size of one frac = 2D (n) + 2D (d) + 2D (e) = 6D.
    #    Total memory z = 2 * 6D = 12D.
    z = 12

    # 6. Print the final answer in the f:z format.
    print(f"{f_rounded}:{z}")

calculate_time_dilation()