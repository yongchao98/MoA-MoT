import math

def solve():
    """
    This script calculates the minimum total heat energy based on the derived analytical solutions.
    """
    
    # The total energy E_total is given by the integral of Phi(T1(sqrt(2)x) + T2(x)) from 0 to 1.
    # From the fractional derivative constraints, we determine Phi(z) = z/2 + 1.
    # So, E_total = integral_0^1 ( (T1(sqrt(2)x) + T2(x))/2 + 1 ) dx
    # E_total = 1 + 1/2 * integral_0^1(T1(sqrt(2)x))dx + 1/2 * integral_0^1(T2(x))dx

    # The integral of the T1 term is calculated analytically.
    # integral_0^1(T1(sqrt(2)x))dx = 7*pi/6 - 10/3
    integral_T1 = 7 * math.pi / 6 - 10 / 3

    # The minimum value of the integral of the T2 term is calculated analytically.
    # min integral_0^1(T2(x))dx = 41/108
    integral_T2_min = 41 / 108
    
    # Now we assemble the final expression for the minimum total energy.
    # E_total = 1 + 1/2 * (7*pi/6 - 10/3) + 1/2 * (41/108)
    # E_total = 1 + 7*pi/12 - 5/3 + 41/216
    # E_total = 7*pi/12 - 2/3 + 41/216
    # E_total = 7*pi/12 - 144/216 + 41/216
    # E_total = 7*pi/12 - 103/216

    E_total = 7 * math.pi / 12 - 103 / 216

    # Print the components of the final equation as requested.
    print("The final equation for the minimum total energy is:")
    print("E_total = 1 + (1/2) * [integral of T1] + (1/2) * [integral of T2]")
    print(f"Value of the constant term: 1")
    print(f"Value of the T1 integral term: 0.5 * {integral_T1:.8f} = {0.5 * integral_T1:.8f}")
    print(f"Value of the T2 integral term: 0.5 * {integral_T2_min:.8f} = {0.5 * integral_T2_min:.8f}")
    print("\nWhich simplifies to:")
    print("E_total = 7*pi/12 - 103/216")
    
    # Print the final numerical value.
    print(f"\nThe minimum total heat energy E_total is: {E_total:.8f}")

solve()
<<<1.35574386>>>