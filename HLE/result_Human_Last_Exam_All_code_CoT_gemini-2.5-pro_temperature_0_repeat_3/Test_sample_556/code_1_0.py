import sys

# This script is designed to be executed in a terminal.
# If you are using a different environment, the print statements will guide you through the logic.

def solve_price():
    """
    This function calculates the price 'p' based on the economic model described.
    """
    # Step 1: Define knowns from the problem statement
    # Probabilities of (effort, output) pairs
    P_e22_y132 = 0.4375
    P_e22_y44 = 0.0625
    P_e10_y60 = 0.0625
    P_e10_y20 = 0.4375

    # Observed effort levels corresponding to high and low signals
    e_H_obs = 22
    e_L_obs = 10

    # Step 2: Deduce the states of the world, s_H and s_L, from y = s * e
    # When e=22, y can be 132 or 44. So s can be 132/22=6 or 44/22=2.
    # When e=10, y can be 60 or 20. So s can be 60/10=6 or 20/10=2.
    s_H = 132 / e_H_obs
    s_L = 44 / e_H_obs
    print("Step 1: Deduce the states of the world from the data.")
    print(f"The high state of the world is s_H = {s_H}")
    print(f"The low state of the world is s_L = {s_L}")
    print("-" * 50)

    # Step 3: Solve for the model's structural parameters (π, λ, γ)
    # The joint probabilities P(e,s) can be expressed as:
    # P(e_H, s_H) = π * P(θ=s_H|s=s_H) = P_e22_y132
    # P(e_H, s_L) = (1-π) * P(θ=s_H|s=s_L) = P_e22_y44
    # P(e_L, s_H) = π * P(θ=s_L|s=s_H) = P_e10_y60
    # P(e_L, s_L) = (1-π) * P(θ=s_L|s=s_L) = P_e10_y20
    # where P(θ|s) are functions of λ and γ.
    # Let A = P(θ=s_H|s=s_H), C = P(θ=s_L|s=s_H). Then A+C=1.
    # From the equations: π*A = 0.4375 and π*C = 0.0625.
    # Dividing them gives A/C = 7. With A+C=1, we get C=1/8 and A=7/8.
    A = 0.4375 / (0.4375 + 0.0625)
    C = 1 - A
    # Let B = P(θ=s_H|s=s_L), D = P(θ=s_L|s=s_L). Then B+D=1.
    # From the equations: (1-π)*B = 0.0625 and (1-π)*D = 0.4375.
    # Dividing them gives B/D = 1/7. With B+D=1, we get B=1/8 and D=7/8.
    B = 0.0625 / (0.0625 + 0.4375)
    D = 1 - B

    # Now, solve for λ and γ using their definitions:
    # B = γ(1-λ) = 1/8
    # C = (1-γ)(1-λ) = 1/8
    # This implies γ = 1-γ, so γ = 0.5.
    # Substituting γ=0.5 into B gives 0.5*(1-λ)=0.125, so 1-λ=0.25, and λ=0.75.
    gamma = 0.5
    lambda_val = 0.75

    # Finally, solve for π using π*A = 0.4375
    pi = P_e22_y132 / A
    
    print("Step 2: Solve for the model's structural parameters.")
    print(f"Probability of high state, π = {pi}")
    print(f"Signal accuracy, λ = {lambda_val}")
    print(f"Independent signal probability, γ = {gamma}")
    print("-" * 50)

    # Step 4: Calculate the employee's conditional expectation of the state E[s|θ]
    # P(θ=s_H) = P(e=e_H) = P(e_H, s_H) + P(e_H, s_L) = 0.4375 + 0.0625 = 0.5
    P_theta_H = P_e22_y132 + P_e22_y44
    # P(s=s_H|θ=s_H) = P(e_H, s_H) / P(θ=s_H)
    P_sH_given_thetaH = P_e22_y132 / P_theta_H
    P_sL_given_thetaH = 1 - P_sH_given_thetaH
    E_s_given_thetaH = s_H * P_sH_given_thetaH + s_L * P_sL_given_thetaH

    # P(θ=s_L) = P(e=e_L) = P(e_L, s_H) + P(e_L, s_L) = 0.0625 + 0.4375 = 0.5
    P_theta_L = P_e10_y60 + P_e10_y20
    # P(s=s_H|θ=s_L) = P(e_L, s_H) / P(θ=s_L)
    P_sH_given_thetaL = P_e10_y60 / P_theta_L
    P_sL_given_thetaL = 1 - P_sH_given_thetaL
    E_s_given_thetaL = s_H * P_sH_given_thetaL + s_L * P_sL_given_thetaL

    print("Step 3: Calculate conditional expectations of the state.")
    print(f"E[s | θ=s_H] = {s_H}*{P_sH_given_thetaH:.3f} + {s_L}*{P_sL_given_thetaH:.3f} = {E_s_given_thetaH}")
    print(f"E[s | θ=s_L] = {s_H}*{P_sH_given_thetaL:.3f} + {s_L}*{P_sL_given_thetaL:.3f} = {E_s_given_thetaL}")
    print("-" * 50)

    # Step 5: Determine β from the incentive compatibility constraint e = β * E[s|θ]
    beta_from_H = e_H_obs / E_s_given_thetaH
    beta_from_L = e_L_obs / E_s_given_thetaL
    # The values must be consistent
    beta = beta_from_H

    print("Step 4: Determine the contract parameter β.")
    print(f"From the high-effort case: β = e_H / E[s|θ=s_H] = {e_H_obs} / {E_s_given_thetaH} = {beta}")
    print(f"From the low-effort case:  β = e_L / E[s|θ=s_L] = {e_L_obs} / {E_s_given_thetaL} = {beta_from_L}")
    print("-" * 50)

    # Step 6: Determine p from the firm's profit maximization problem.
    # The firm's optimization problem implies p = β.
    p = beta

    print("Step 5: Determine the price p.")
    print("The firm's profit maximization implies a simple relationship between p and β.")
    print("Final Equation: p = β")
    print(f"p = {beta}")
    print("-" * 50)
    
    return p

if __name__ == '__main__':
    price = solve_price()
    # The final answer is printed here, enclosed in <<<>>>
    # To run this script, save it as a .py file and execute from the terminal.
    # For example: python your_script_name.py
    # The script will print the step-by-step derivation and then the final answer.
    final_answer = int(price) if price.is_integer() else price
    sys.stdout.write(f"<<<{final_answer}>>>")