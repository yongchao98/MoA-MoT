import math

def solve_random_walk_probability():
    """
    Calculates the probability that a 2D random walk conditioned to avoid the origin,
    starting from (3000, 4000), will never hit the set of the origin's four neighbors.
    """
    # Starting point coordinates
    x_coord = 3000
    y_coord = 4000

    # Calculate the distance r = |x_0| from the origin
    r = math.sqrt(x_coord**2 + y_coord**2)

    # The potential kernel at a neighbor of the origin, a(1,0), has an exact value.
    # a(1,0) = 4/pi - 1
    a1 = 4 / math.pi - 1

    # For a large distance r, the potential kernel a(r) is given by the asymptotic formula:
    # a(r) ~ (2/pi) * log(r) + (2*gamma/pi) + (1/pi) * log(8)
    # where gamma is the Euler-Mascheroni constant.
    gamma = math.eulgamma
    a_r = (1 / math.pi) * (2 * math.log(r) + 2 * gamma + math.log(8))

    # The probability of the conditioned walk hitting the neighbors of the origin is P_hit ≈ a(1,0) / a(x_0)
    p_hit = a1 / a_r

    # The probability of never hitting the set is 1 - P_hit
    p_avoid = 1 - p_hit

    print("Problem: Find the probability that a 2D random walk conditioned to avoid the origin,")
    print(f"starting from ({x_coord},{y_coord}), never hits the origin's four neighbors.")
    print("\n--- Calculation Steps ---")
    print(f"1. Starting distance from origin: r = sqrt({x_coord}^2 + {y_coord}^2) = {r}")
    print(f"2. Potential kernel at a neighbor of the origin: a(1,0) = 4/π - 1 ≈ {a1:.5f}")
    print(f"3. Potential kernel at the starting point (asymptotic formula): a({r:.0f}) ≈ {a_r:.5f}")
    
    print("\n--- Final Equation ---")
    print("The probability of avoiding the neighbors is P_avoid ≈ 1 - a(1,0) / a(r)")
    print(f"P_avoid ≈ 1 - {a1:.5f} / {a_r:.5f}")
    print(f"P_avoid ≈ 1 - {p_hit:.5f}")
    print(f"P_avoid ≈ {p_avoid:.5f}")

    # Format the final answer to two significant digits.
    # For a number like 0.957..., the first two significant digits are 9 and 5.
    # The third digit (7) causes rounding up, so the result is 0.96.
    final_answer_formatted = f"{p_avoid:.2g}"
    print(f"\nThe approximate answer with two significant digits is: {final_answer_formatted}")

solve_random_walk_probability()
<<<0.96>>>