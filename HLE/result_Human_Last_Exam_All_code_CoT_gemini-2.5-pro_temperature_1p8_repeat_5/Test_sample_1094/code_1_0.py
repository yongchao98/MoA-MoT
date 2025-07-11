import numpy as np
from scipy.optimize import brentq

def calculate_ac_loss(i, aspect_ratio, verbose=True):
    """
    Calculates the normalized AC loss 2*pi*Q/(mu_0*Ic^2) for a superconducting
    elliptic bar with an AC transport current.

    The calculation is based on the model by Clem and Sanchez, JAP 114, 023903 (2013).

    Args:
        i (float): Normalized current amplitude, i = Im/Ic. Must be 0 < i < 1.
        aspect_ratio (float): The aspect ratio k = b/a of the ellipse's semi-axes.
                              Must be 0 < k <= 1.
        verbose (bool): If True, prints the intermediate steps of the calculation.

    Returns:
        float: The normalized AC loss.
    """
    if not (0 < i < 1):
        raise ValueError("Current ratio i must be between 0 and 1 (exclusive).")
    if not (0 < aspect_ratio <= 1):
        raise ValueError("Aspect ratio b/a must be between 0 (exclusive) and 1 (inclusive).")

    k = aspect_ratio

    if np.isclose(k, 1.0):
        # Special case: round wire (k=1)
        # Using the well-established formula for bipolar current loss in a round wire
        if verbose:
            print(f"Calculating for round wire (k=1.0) with i = {i}")
        if i < 0.5:
            # For i < 0.5, loss from -Im to Im is same as 0 to 2Im
            arg = 2*i
            loss_val = 2 * (arg - arg**2/2 - (1-arg)*np.log(1-arg))
            norm_loss = 2 * loss_val
        else: # 0.5 <= i < 1
            # For i >= 0.5, the loss cycle is more complex
            # Q_N_bipolar = Q_N_mono(1) + Q_N_mono(2i-1)
            q_n_mono_1 = 0.5
            arg = 2*i-1
            q_n_mono_2i_1 = arg - arg**2/2 - (2-2*i)*np.log(2-2*i)
            q_n_bipolar = q_n_mono_1 + q_n_mono_2i_1
            norm_loss = 2 * q_n_bipolar
        if verbose:
             print(f"Final Normalized Loss = {norm_loss:.6f}")
        return norm_loss


    # General case for an ellipse (k < 1)
    e2 = 1.0 - k**2
    e = np.sqrt(e2)

    if verbose:
        print(f"Calculating for ellipse with i = {i}, aspect ratio k = b/a = {k}")
        print(f"Eccentricity e = sqrt(1-k^2) = {e:.6f}")

    # --- Find penetration parameters g_m and g_p ---
    # These parameters characterize the current-free region inside the conductor.
    
    # 1. Solve for g_m from i = 1 - (g_m/k) * sqrt(g_m^2 - e^2)
    def f_gm(g):
        # We need to find the root of this function for g in [e, 1]
        return (g / k) * np.sqrt(g**2 - e2) - (1.0 - i)
    g_m = brentq(f_gm, e, 1.0)

    # 2. Solve for g_p from (g_p/k) * sqrt(g_p^2 - e^2) = 1 + i
    def f_gp(g):
        # We need to find the root of this function for g > e
        return (g / k) * np.sqrt(g**2 - e2) - (1.0 + i)
    # The upper bound 10 is arbitrary but safely larger than the expected root.
    g_p = brentq(f_gp, e, 10.0)

    if verbose:
        print(f"Solved for penetration parameters: g_m = {g_m:.6f}, g_p = {g_p:.6f}")

    # --- Calculate loss using the general formula L = L_strip - L_correction ---
    
    # 3. Loss for a thin strip (limit k->0)
    L_strip = 2.0 * ((1+i)*np.log(1+i) + (1-i)*np.log(1-i) - i**2)

    # 4. Correction term F_e for the elliptical shape
    # This term corrects the thin strip formula for finite aspect ratio
    term1 = (2.0/g_m) * np.arctanh(e/g_m)
    term2 = (2.0/g_p) * np.arctanh(e/g_p)

    log_term1_num = g_m**2 - e2
    log_term1_den = g_p**2 - e2
    log_term1 = np.log(log_term1_num / log_term1_den)
    
    log_term2_num = g_m**2 * (1 - e2 * g_p**2)
    log_term2_den = g_p**2 * (1 - e2 * g_m**2)
    log_term2 = np.log(log_term2_num / log_term2_den)

    factor1 = (1 - i - k**2/g_m**2)
    
    F_e = term1 - term2 + factor1 * log_term1 + 2*i * log_term2
    
    # 5. Full correction term L_corr
    L_corr = (2.0 * k**2 / e**4) * F_e

    # 6. Final normalized loss
    final_loss = L_strip - L_corr

    if verbose:
        print("\n--- Final Equation ---")
        print("Loss = L_strip - L_correction")
        print(f"L_strip({i:.2f}) = {L_strip:.6f}")
        print(f"L_correction = (2*k^2/e^4) * F_e = (2*{k**2:.2f}/{e**4:.4f}) * {F_e:.6f} = {L_corr:.6f}")
        print(f"Final Normalized Loss = {L_strip:.6f} - {L_corr:.6f} = {final_loss:.6f}")

    return final_loss

if __name__ == '__main__':
    # --- Example Calculation ---
    # User can modify these values
    current_ratio_i = 0.6
    aspect_ratio_k = 0.5  # This is b/a

    # Calculate and print the result with intermediate steps
    normalized_loss = calculate_ac_loss(current_ratio_i, aspect_ratio_k)
    
    # The final answer in the requested format
    # Note: We provide the script and its output. The final numerical value is part of the output.
    print(f"\n<<<The normalized AC loss is {normalized_loss}>>>")
