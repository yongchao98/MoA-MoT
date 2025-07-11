import math

def solve():
    """
    Calculates the minimum q for Alice's winning probability to be at least 95%.
    """
    # The desired probability for Alice to win. This is 'x' in our derivation.
    p_alice_win = 0.95

    # From the derived formula: q^4 * x^3 = 1 - (1-x)^(1/3)
    # We solve for q.

    # 1. Calculate the right-hand side (RHS) of the equation
    # (1-x) is the probability that the game does NOT result in a win for Alice.
    # In our problem, this corresponds to the value 0.05.
    one_minus_x = 1 - p_alice_win
    rhs = 1 - (one_minus_x)**(1/3)

    # 2. Calculate the x^3 term on the left-hand side (LHS)
    x_cubed = p_alice_win**3

    # 3. Solve for q^4
    q0_pow4 = rhs / x_cubed

    # 4. Solve for q0, the minimum value of q
    q0 = q0_pow4**(1/4)

    # As requested, printing each number in the final equation.
    print("Derived relationship between Alice's win probability (x) and edge probability (q):")
    print("q^4 * x^3 = 1 - (1 - x)^(1/3)\n")
    print(f"Setting x = {p_alice_win}:")
    print(f"q0^4 * ({p_alice_win})^3 = 1 - (1 - {p_alice_win})^(1/3)")
    print(f"q0^4 * {x_cubed:.6f} = 1 - ({one_minus_x})^(1/3)")
    print(f"q0^4 * {x_cubed:.6f} = 1 - {one_minus_x**(1/3):.6f}")
    print(f"q0^4 * {x_cubed:.6f} = {rhs:.6f}")
    print(f"q0^4 = {rhs:.6f} / {x_cubed:.6f}")
    print(f"q0^4 = {q0_pow4:.6f}")
    print(f"q0 = ({q0_pow4:.6f})^(1/4) = {q0:.6f}\n")

    # 5. Calculate the final answer required by the problem: floor(100 * q0)
    final_answer = math.floor(100 * q0)

    print(f"The minimum value q0 is approximately {q0:.6f}.")
    print(f"The value of 100 * q0 is {100 * q0:.6f}.")
    print(f"The final answer, floor(100 * q0), is: {final_answer}")


solve()