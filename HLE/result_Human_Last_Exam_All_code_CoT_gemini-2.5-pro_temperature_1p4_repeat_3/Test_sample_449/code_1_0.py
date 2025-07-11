import math

def solve():
    """
    Calculates the probability that the conditioned random walk never enters the set 
    of the four neighbors of the origin.
    """
    
    # Constants
    PI = math.pi
    EULER_GAMMA = 0.5772156649  # Euler-Mascheroni constant

    # Starting position and its magnitude
    x0 = (3000, 4000)
    r = math.sqrt(x0[0]**2 + x0[1]**2)

    # Calculate a_1 = a(1,0), the potential kernel at a neighbor of the origin.
    # The exact value is a_1 = 8/pi - 5/3.
    a1_val = 8 / PI - 5 / 3
    print(f"The potential kernel at a neighbor of the origin is a_1 = 8/\u03C0 - 5/3 \u2248 {a1_val:.4f}")

    # Calculate the constant C in the asymptotic expansion of a(x).
    # C = (2*gamma + ln(8))/pi
    C = (2 * EULER_GAMMA + math.log(8)) / PI
    print(f"The constant in the asymptotic expansion is C = (2\u03B3 + ln(8))/\u03C0 \u2248 {C:.4f}")

    # Calculate a(x_0) using the asymptotic formula for large |x_0|.
    # a(x_0) ~ (2/pi)*ln|x_0| + C
    ax0_val = (2 / PI) * math.log(r) + C
    print(f"The potential kernel at x_0=(3000,4000) is a(x_0) \u2248 (2/\u03C0)*ln({r:.0f}) + C \u2248 {ax0_val:.4f}")

    # Calculate the final probability p(x_0) = 1 - a_1/a(x_0).
    prob = 1 - a1_val / ax0_val
    
    print("\nThe probability is given by the formula p = 1 - a_1/a(x_0).")
    print("Plugging in the computed values:")
    print(f"p = 1 - {a1_val:.4f} / {ax0_val:.4f} = {prob:.4f}")

    # Final answer with two significant digits
    final_answer = round(prob, 2)
    print(f"\nThe approximate answer with two significant digits is: {final_answer}")

solve()