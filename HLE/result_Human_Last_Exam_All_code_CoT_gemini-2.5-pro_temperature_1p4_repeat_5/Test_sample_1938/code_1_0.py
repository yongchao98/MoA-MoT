import math

def solve():
    """
    This function calculates the minimum q for Alice to win with at least 95% probability.
    """
    # Step 1: Define the relationship between Alice's winning probability p and the edge probability q.
    # The derived equation is p = 1 - (1 - q^4 * p^3)^3.
    # For a detailed derivation, please refer to the explanation above.

    # Step 2: Set the desired winning probability for Alice.
    p = 0.95
    print(f"We are looking for the minimum q such that Alice's winning probability, p, is at least 0.95.")
    print(f"We solve for q when p = {p}.")
    print(f"The equation relating p and q is: p = 1 - (1 - q^4 * p^3)^3")
    print("-" * 30)

    # Step 3: Substitute p = 0.95 into the equation and solve for q.
    print(f"Substitute p = {p} into the equation:")
    print(f"{p} = 1 - (1 - q^4 * ({p})^3)^3")

    # Isolate the term with q
    p_cubed = p**3
    print(f"\nRearranging the equation to solve for q:")
    print(f"(1 - q^4 * {p_cubed:.6f})^3 = 1 - {p}")
    
    val_1_minus_p = 1 - p
    print(f"(1 - q^4 * {p_cubed:.6f})^3 = {val_1_minus_p}")

    # Take the cube root of both sides
    val_cbrt = val_1_minus_p**(1/3)
    print(f"\nTake the cube root of both sides:")
    print(f"1 - q^4 * {p_cubed:.6f} = ({val_1_minus_p})^(1/3) = {val_cbrt:.6f}")

    # Isolate the q^4 term
    val_numerator = 1 - val_cbrt
    print(f"\nContinue rearranging:")
    print(f"q^4 * {p_cubed:.6f} = 1 - {val_cbrt:.6f} = {val_numerator:.6f}")

    # Solve for q^4
    q_fourth = val_numerator / p_cubed
    print(f"\nSolve for q^4:")
    print(f"q^4 = {val_numerator:.6f} / {p_cubed:.6f} = {q_fourth:.6f}")

    # Solve for q
    q0 = q_fourth**(1/4)
    print(f"\nSolve for q by taking the fourth root:")
    print(f"q = ({q_fourth:.6f})^(1/4) = {q0:.6f}")
    print("-" * 30)

    # Step 4: Calculate the final answer as requested.
    final_answer = math.floor(100 * q0)
    print(f"The minimum value of q is q0 ≈ {q0:.6f}")
    print(f"The question asks for the value of floor(100 * q0).")
    print(f"100 * q0 ≈ {100 * q0:.6f}")
    print(f"The final result is floor({100 * q0:.6f}) = {final_answer}")

if __name__ == '__main__':
    solve()