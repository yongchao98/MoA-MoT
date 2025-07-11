import math

def solve_contracts():
    """
    This function calculates the separating equilibrium contracts for Chelsea and the handymen.
    """
    # Step 1: Solve for the Low-Type's First-Best Contract (v_L, p_L)
    # Chelsea maximizes v_L - p_L subject to p_L - v_L^2 >= 0.
    # This is equivalent to maximizing v_L - v_L^2.
    # The derivative is 1 - 2*v_L. Setting to zero gives the optimal v_L.
    # 1 - 2*v_L = 0
    v_L = 0.5
    
    # The productivity constraint is v_L <= 20 - (1/3)*v_L => (4/3)*v_L <= 20 => v_L <= 15.
    # Since 0.5 <= 15, the value is valid.
    
    # The price p_L is set to make the low-type's participation constraint bind (u_L = 0).
    # p_L - v_L^2 = 0 => p_L = v_L^2
    p_L = v_L**2

    # Step 2: Solve for the High-Type's First-Best Contract (v_H, p_H)
    # Chelsea maximizes v_H - p_H subject to p_H - v_H >= 0.
    # She sets p_H = v_H, making her profit v_H - p_H = 0.
    # As a tie-breaker, she chooses the highest possible v_H to maximize her house value.
    # The productivity constraint is v_H <= 30.
    v_H = 30.0
    
    # The price p_H is set to make the high-type's participation constraint bind (u_H = 0).
    # p_H - v_H = 0 => p_H = v_H
    p_H = v_H

    # Step 3: Output the results
    # The final answer is the pair of contracts.
    print("The separating equilibrium consists of the following pair of contracts:")
    print(f"Contract for the Low-Type Handyman: (v_L, p_L) = ({v_L}, {p_L})")
    print(f"Contract for the High-Type Handyman: (v_H, p_H) = ({v_H}, {p_H})")

solve_contracts()