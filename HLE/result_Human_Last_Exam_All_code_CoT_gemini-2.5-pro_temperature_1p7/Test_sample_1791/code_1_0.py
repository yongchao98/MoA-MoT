def get_effective_interaction():
    """
    This function derives and prints the formula for the effective 
    electron-electron interaction mediated by phonons.
    """
    
    # Define symbolic variables for the formula as strings
    g = "g"          # electron-phonon coupling constant
    q_squared = "|q|^2"  # squared magnitude of the phonon wave vector
    m = "m"          # ionic mass
    w_q = "w_q"      # phonon frequency for wave vector q
    omega = "w"      # electron energy transfer
    n_q = "n_q"      # electron density operator for wave vector q
    n_mq = "n_{-q}"  # electron density operator for wave vector -q

    # The effective interaction for a single mode q is V_eff(q) * n_q * n_{-q}
    # (The factor of 1/2 appears when summing over all q to avoid double counting)
    # Here, we present the term for a specific q.
    
    print("The effective electron-electron interaction term mediated by phonons for a mode 'q' is:")
    
    # We construct the formula string part by part as requested.
    # The potential is V_eff(q, w) = (g^2 * |q|^2) / (m * (w_q^2 - w^2))
    # It corresponds to an attractive interaction when w < w_q.
    
    # Numerator of the potential V_eff
    numerator_V = f"{g}^2 * {q_squared}"
    
    # Denominator of the potential V_eff
    denominator_V = f"{m} * ({w_q}^2 - {omega}^2)"

    # Print the full interaction term: V_eff * n_q * n_{-q}
    # Note: A minus sign is conventionally included to show attraction explicitly, 
    # making the potential -g^2|q|^2/(m(w^2-w_q^2)). We use the derived form.
    print(f"({numerator_V}) / ({denominator_V}) * {n_q} * {n_mq}")

    print("\nWhere:")
    print(f"  g: electron-phonon coupling constant")
    print(f"  q: phonon wave vector")
    print(f"  m: ionic mass")
    print(f"  w_q: phonon frequency")
    print(f"  w: electron energy transfer")
    print(f"  n_q: electron density operator for wavevector q")


if __name__ == '__main__':
    get_effective_interaction()