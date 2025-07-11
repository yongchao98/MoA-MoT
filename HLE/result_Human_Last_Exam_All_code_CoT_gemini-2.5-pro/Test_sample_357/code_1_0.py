import math

def solve_contract_problem():
    """
    This function calculates the separating equilibrium contracts based on the problem's parameters.
    """
    # 1. Define probabilities
    prob_h = 5/6
    prob_l = 1 - prob_h

    # 2. Determine the optimal value for the low type, v_L.
    # Chelsea's profit simplifies to (1/6) * (v_L - v_L^2).
    # To maximize this, we find the vertex of the parabola f(v) = v - v^2.
    # The derivative is 1 - 2*v. Setting to 0 gives v = 0.5.
    opt_v_l = 0.5

    # 3. Determine the optimal value for the high type, v_H.
    # The profit function is independent of v_H. We apply the "no distortion at the top"
    # principle and set v_H to its maximum possible value, given by its productivity.
    theta_h = 30
    opt_v_h = float(theta_h)

    # 4. Calculate the prices based on the binding constraints.
    # For this optimal case, the Individual Rationality (IR) constraints for both types are binding.
    # IR-L: p_L - v_L^2 = 0  => p_L = v_L^2
    # IR-H: p_H - v_H = 0    => p_H = v_H
    opt_p_l = opt_v_l**2
    opt_p_h = opt_v_h

    # 5. Print the results clearly.
    print("Based on the profit maximization for Chelsea subject to the handymen's constraints, the separating equilibrium is found.")
    print("The contract for the low type (L) is:")
    print(f"(v_L, p_L) = ({opt_v_l}, {opt_p_l})")
    print("\nThe contract for the high type (H) is:")
    print(f"(v_H, p_H) = ({opt_v_h}, {opt_p_h})")

    # This part is for the final answer format and is not printed to the user.
    global final_answer
    final_answer = f"<<<(({opt_v_l}, {opt_p_l}), ({opt_v_h}, {opt_p_h}))>>>"

# Run the solver
solve_contract_problem()