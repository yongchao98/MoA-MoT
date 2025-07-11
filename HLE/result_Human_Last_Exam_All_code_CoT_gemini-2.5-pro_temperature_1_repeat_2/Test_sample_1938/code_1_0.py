import math

def solve():
    """
    This function solves for the minimum q such that Alice's win probability is >= 95%.
    """
    
    # We need to solve the cubic equation for z = q^4 derived from the game's recurrence relations.
    # The equation is p^2*z^3 - 3*p*z^2 + 3*z - 1 = 0, where p is Alice's win probability.
    # For p = 0.95, the equation becomes:
    # 0.9025 * z^3 - 2.85 * z^2 + 3*z - 1 = 0.
    # The numbers in this final equation are:
    a = 0.9025  # p^2
    b = -2.85   # -3*p
    c = 3.0     # 3
    d = -1.0    # -1

    def g(z):
        """The cubic function we need to find the root of."""
        return a * z**3 + b * z**2 + c * z + d

    # We use binary search to find the root z.
    # A non-trivial solution exists for q^4 > 1/3.
    low_z = 1.0 / 3.0
    high_z = 1.0
    
    # Iterate 100 times for high precision.
    for _ in range(100):
        mid_z = (low_z + high_z) / 2
        if g(mid_z) < 0:
            low_z = mid_z
        else:
            high_z = mid_z
    
    z_0 = (low_z + high_z) / 2
    
    # Now, calculate q_0 from z_0 = q_0^4
    q_0 = z_0**0.25
    
    # The final answer is floor(100 * q_0)
    result = math.floor(100 * q_0)
    
    print(result)

solve()