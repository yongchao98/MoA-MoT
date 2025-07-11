import numpy as np
from scipy.special import eval_genlaguerre, gammaln

def calculate_ratio_components(n):
    """
    Calculates the components of the ratio D_n(r*) / D_n^c(r*) for a given n.
    All calculations are done in atomic units (a_0 = 1, e = 1, m_e = 1, hbar = 1).
    """
    if n <= 0:
        return None

    # 1. Classical part
    # r* is the radius that maximizes the classical distribution D_n^c(r)
    r_star = 1.5 * n**2
    # Value of the normalized classical distribution at r*
    D_c_peak_val = 9 / (2 * np.pi * np.sqrt(3) * n**2)

    # 2. Quantum part
    # We need to calculate D_n(r*) = (r*)^2 * sum_{l=0}^{n-1} (2l+1) * |R_nl(r*)|^2
    quantum_sum = 0
    rho = 2 * r_star / n

    for l in range(n):
        # R_nl^2 = N_nl^2 * exp(-rho) * rho^(2l) * (L_{n-l-1}^{2l+1}(rho))^2
        # N_nl^2 = (2/n)^3 * (n-l-1)! / (2n * (n+l)!)
        # Use log-gamma for numerical stability of factorials: log(k!) = gammaln(k+1)
        log_N_nl_sq = 3 * np.log(2/n) + gammaln(n-l) - np.log(2*n) - gammaln(n+l+1)
        N_nl_sq = np.exp(log_N_nl_sq)
        
        # Associated Laguerre polynomial L_{n-l-1}^{2l+1}(rho)
        lag_val = eval_genlaguerre(n - l - 1, 2 * l + 1, rho)

        R_nl_sq = N_nl_sq * np.exp(-rho) * (rho**(2*l)) * (lag_val**2)
        quantum_sum += (2*l + 1) * R_nl_sq
        
    D_q_val = r_star**2 * quantum_sum
    
    # 3. The ratio
    ratio = D_q_val / D_c_peak_val
    
    return {
        "n": n,
        "r_star": r_star,
        "D_q_val_expr": f"{D_q_val/np.exp(-r_star/n):.4f} * exp(-{r_star/n})",
        "D_q_val": D_q_val,
        "D_c_peak_val_expr": f"9 / (2 * pi * sqrt(3) * {n**2})",
        "D_c_peak_val": D_c_peak_val,
        "ratio": ratio
    }

def solve():
    """
    Finds the maximum ratio by calculating it for n from 1 to 20.
    """
    max_ratio_info = {"ratio": 0}
    
    for n in range(1, 21):
        current_info = calculate_ratio_components(n)
        if current_info["ratio"] > max_ratio_info["ratio"]:
            max_ratio_info = current_info

    n_max = max_ratio_info['n']
    r_star_max = max_ratio_info['r_star']
    D_q_max = max_ratio_info['D_q_val']
    D_c_max = max_ratio_info['D_c_peak_val']
    ratio_max = max_ratio_info['ratio']
    
    # Exact symbolic expression for the ratio at n=3
    # D_q(n=3) = (25581/4) * exp(-9)
    # D_c(n=3) = 1 / (2*pi*sqrt(3))
    # Ratio = D_q/D_c = (25581/4) * exp(-9) * 2 * pi * sqrt(3) = 12790.5 * pi * sqrt(3) * exp(-9)

    print(f"The maximum ratio occurs at n = {n_max}.")
    print("\n--- Details of the calculation for n = {} ---".format(n_max))
    print(f"The classical distribution is maximized at r* = {r_star_max:.1f} a_0.")
    print(f"The quantum distribution at this radius is D_n(r*) = {D_q_max:.8e}.")
    print(f"The classical distribution at this radius is D_n^c(r*) = {D_c_max:.8e}.")
    print("The final equation for the maximum ratio is:")
    print(f"Ratio = D_{n_max}({r_star_max:.1f}) / D_{n_max}^c({r_star_max:.1f}) = {D_q_max:.4e} / {D_c_max:.4e}")
    print(f"The value of the maximum ratio is {ratio_max:.8f}.")

solve()
<<<8.60450091>>>