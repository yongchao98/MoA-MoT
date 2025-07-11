import math

def solve_integral():
    """
    This function calculates the number of poles inside the contour Gamma,
    and then computes the value of the contour integral.
    """
    num_poles = 0
    pi = math.pi

    # From theoretical analysis of the curve equation, we can deduce that
    # the curve is bounded. The relevant integer values for k are in a small
    # range around 0, and for a are in a small range around 0. We'll use a
    # search space that safely covers this region.
    k_range = range(-4, 5)
    a_range = range(-40, 41)

    for k in k_range:
        for a in a_range:
            p = float(a + 2 * pi * k)
            q = float(a - 2 * pi * k)

            # Necessary condition for a point to be on/inside the curve
            # is (p^2) <= 573.77. We use 574 as a safe upper bound.
            if p**2 > 574:
                # Can break early if p is monotonic in this 'a' loop
                if (a > 0 and 2 * pi * k >= 0) or (a < 0 and 2*pi*k <= 0):
                  if p**2 > (a_range.start)**2 : # if p is moving away from 0
                    break

            # The pole z = a + 2*pi*i*k is inside Gamma if the following expression is < 0.
            # This is the simplified equation for the curve Gamma in (p,q) coordinates.
            # 1200*q^4 - 20*p^3*q^2 + 3*(p^2 - 400)^3 < 0
            
            p2 = p * p
            p3 = p2 * p
            q2 = q * q
            q4 = q2 * q2

            value = 1200 * q4 - 20 * p3 * q2 + 3 * (p2 - 400)**3

            if value < 0:
                num_poles += 1
    
    # The value of the integral is 2 * pi * i * (number of poles inside Gamma).
    # The residue at each pole is 1.
    
    coeff = 2 * num_poles
    
    print(f"The number of enclosed poles is N_0 = {num_poles}.")
    print(f"The integral is 2 * pi * i * N_0.")
    print(f"The equation for the integral is: 2 * pi * i * {num_poles} = {coeff} * pi * i")

    # The final answer in the requested format
    print("\n<<<{}*pi*i>>>".format(coeff))

solve_integral()