import sys

def solve_mass_problem(n, q, q_v):
    """
    Calculates the total mass based on the provided parameters.

    Args:
      n (int): The dimension for the general linear group.
      q (int): The characteristic of the local field (size of the constant field).
      q_v (int): The order of the residual field.
    """
    # The problem asks for the total mass of (q_v * (q - 1)) / (q_v - 1) * mu.
    # Our analysis shows that the mass of the measure mu on the space of lattices is 1 / (q - 1).
    # This mass seems to be independent of n.

    if not isinstance(q, int) or q <= 1:
        print("Error: q must be an integer greater than 1.", file=sys.stderr)
        return
    if not isinstance(q_v, int) or q_v <= 1:
        print("Error: q_v must be an integer greater than 1.", file=sys.stderr)
        return

    # mu_mass = 1 / (q - 1)
    # The expression to calculate is (q_v * (q - 1) / (q_v - 1)) * mu_mass
    # This simplifies to q_v / (q_v - 1)

    result = q_v / (q_v - 1)

    print("The total mass is given by the expression:")
    print(f"  (q_v * (q - 1)) / (q_v - 1) * mu")
    print("\nBased on our analysis, the mass mu is independent of n and equals 1 / (q - 1).")
    print("Substituting this value:")
    print(f"  Mass = ({q_v} * ({q} - 1)) / ({q_v} - 1) * (1 / ({q} - 1))")
    print("\nThe terms involving (q - 1) cancel out, simplifying the expression to:")
    print(f"  Mass = {q_v} / ({q_v} - 1)")
    print("\nFinally, computing the numerical value:")
    print(f"  Mass = {q_v} / {q_v - 1} = {result}")

# Example values for n, q, and q_v.
# You can change these values to see the result for a different case.
n_val = 2
q_val = 5
q_v_val = 25

solve_mass_problem(n_val, q_val, q_v_val)
