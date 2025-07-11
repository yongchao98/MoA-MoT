import math

def solve():
    """
    Calculates the minimum q for Alice to win with at least 95% probability.
    """
    # The problem asks for the minimum q such that Alice's winning probability, p_A, is >= 0.95.
    # We solve for the boundary condition where p_A is exactly 0.95.
    p_A = 0.95

    print(f"We are looking for the minimum q where Alice's win probability, p_A, is at least 0.95.")
    print(f"We solve the system of equations for the boundary case p_A = {p_A}.")
    print("-" * 30)
    print("The equations are:")
    print("1) p_A = 1 - (1 - q * p_B)^3")
    print("2) p_B = (q * p_A)^3")
    print("-" * 30)

    # Step 1: From equation (1), express q * p_B in terms of p_A.
    # 1 - p_A = (1 - q * p_B)^3
    # (1 - p_A)^(1/3) = 1 - q * p_B
    # q * p_B = 1 - (1 - p_A)^(1/3)
    
    # For clarity in calculation, let's define an intermediate variable y.
    y = (1 - p_A)**(1/3)
    
    print(f"Let's define y = (1 - p_A)^(1/3).")
    print(f"y = (1 - {p_A})^(1/3) = {0.05}^(1/3) = {y:.6f}")
    print("\nThen, q * p_B = 1 - y")
    print("-" * 30)

    # Step 2: Find p_B in terms of p_A.
    # From p_B = q^3 * p_A^3, we have p_B^(1/3) = q * p_A.
    # Dividing the two expressions: (q * p_B) / (q * p_A) = (1 - y) / p_B^(1/3)
    # p_B / p_A = (1 - y) / p_B^(1/3)
    # p_B^(4/3) = p_A * (1 - y)
    # p_B = (p_A * (1 - y))^(3/4)
    p_B = (p_A * (1 - y))**(3/4)

    print("Now we solve for p_B in terms of p_A:")
    print("p_B = (p_A * (1 - y))^(3/4)")
    print(f"p_B = ({p_A} * (1 - {y:.6f}))^(3/4)")
    print(f"p_B = ({p_A * (1 - y):.6f})^(3/4) = {p_B:.6f}")
    print("-" * 30)

    # Step 3: Solve for q_0, the value of q when p_A = 0.95.
    # q_0 = (1 - y) / p_B
    q_0 = (1 - y) / p_B

    print("Finally, we solve for q_0:")
    print("q_0 = (1 - y) / p_B")
    print(f"q_0 = (1 - {y:.6f}) / {p_B:.6f}")
    print(f"q_0 = {1 - y:.6f} / {p_B:.6f} = {q_0:.6f}")
    print("-" * 30)
    
    # Step 4: Calculate the final answer as floor(100 * q_0).
    final_answer = math.floor(100 * q_0)

    print("The problem asks for floor(100 * q_0):")
    print(f"100 * q_0 = 100 * {q_0:.6f} = {100 * q_0:.6f}")
    print(f"floor(100 * q_0) = floor({100 * q_0:.6f}) = {final_answer}")
    
    return final_answer

if __name__ == "__main__":
    result = solve()
    print("\n<<<" + str(result) + ">>>")