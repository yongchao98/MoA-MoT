import math

def F_func(u, u_bar):
    """Calculates the flux F."""
    # Ensure u is within [0,1] to avoid domain errors if any
    u = max(0.0, min(1.0, u))
    return u * (1 - u)**2 * math.exp(-u_bar)

def C1_func(u):
    """Calculates the coefficient (4 - 16u + 24u^2 - 12u^3).
    This arises from the derivation."""
    return 4 - 16*u + 24*u**2 - 12*u**3

def C2_func(u):
    """Calculates the coefficient (4 - 6u), related to d(F_11)/d(u_bar)."""
    return 4 - 6*u

def main():
    """
    Calculates the maximum of the target expression by evaluating it
    at the optimal point found through mathematical analysis.
    """
    # The optimal values are found by analyzing the expression for Q.
    # The maximum is achieved with a step-function-like profile for u(x),
    # where u=0 for x<x_0+1 and u=1 for x>=x_0+1.
    # Evaluating at x=x_0 gives:
    u_opt = 0.0
    u_sh_opt = 1.0  # u_sh = u(x+1)
    u_bar_opt = 0.0 # u_bar = integral from x to x+1 of u=0
    u_bar_sh_opt = 1.0 # u_bar_sh = integral from x+1 to x+2 of u=1

    # The simplified expression for the quantity to maximize is:
    # Q = C1(u) * exp(-2*u_bar) * (u_sh - u) + C2(u) * exp(-u_bar) * (F(u, u_bar) - F(u_sh, u_bar_sh))

    # Calculate all intermediate components for clarity
    c1_val = C1_func(u_opt)
    c2_val = C2_func(u_opt)
    f_val = F_func(u_opt, u_bar_opt)
    f_sh_val = F_func(u_sh_opt, u_bar_sh_opt)
    
    exp_term1_factor = math.exp(-2 * u_bar_opt)
    diff_u = u_sh_opt - u_opt
    term1 = c1_val * exp_term1_factor * diff_u

    exp_term2_factor = math.exp(-u_bar_opt)
    diff_f = f_val - f_sh_val
    term2 = c2_val * exp_term2_factor * diff_f
    
    max_q_value = term1 + term2

    # Output the logic and each number in the final equation
    print("The quantity to maximize is Q = (d/dt + F1 * d/dx) * F11.")
    print("After simplification, this becomes:")
    print("Q = (4 - 16u + 24u^2 - 12u^3) * e^(-2*u_bar) * (u_sh - u) + (4 - 6u) * e^(-u_bar) * (F(u, u_bar) - F(u_sh, u_bar_sh))\n")
    print(f"The maximum is achieved for the state (u, u_sh, u_bar, u_bar_sh) = ({u_opt}, {u_sh_opt}, {u_bar_opt}, {u_bar_sh_opt}).\n")
    
    print("Evaluating the first term: (4 - 16u + 24u^2 - 12u^3) * e^(-2*u_bar) * (u_sh - u)")
    print(f"  = ({C1_func(u_opt)}) * e^(-2*{u_bar_opt}) * ({u_sh_opt} - {u_opt})")
    print(f"  = {c1_val} * {exp_term1_factor} * {diff_u} = {term1}\n")

    print("Evaluating the second term: (4 - 6u) * e^(-u_bar) * (F(u, u_bar) - F(u_sh, u_bar_sh))")
    print(f"  F(u={u_opt}, u_bar={u_bar_opt}) = {u_opt}*(1-{u_opt})^2*e^(-{u_bar_opt}) = {f_val}")
    print(f"  F(u_sh={u_sh_opt}, u_bar_sh={u_bar_sh_opt}) = {u_sh_opt}*(1-{u_sh_opt})^2*e^(-{u_bar_sh_opt}) = {f_sh_val}")
    print(f"  = ({C2_func(u_opt)}) * e^(-{u_bar_opt}) * ({f_val} - {f_sh_val})")
    print(f"  = {c2_val} * {exp_term2_factor} * {diff_f} = {term2}\n")

    print("Final Result:")
    print(f"Max Value = {term1} + {term2} = {max_q_value}")

if __name__ == "__main__":
    main()