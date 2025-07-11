import math

def display_bound_formula():
    """
    This function prints the formula for the upper-bound factor of ||B Q_{0, M}||_oo
    and explains its components.
    """
    
    # The components of the formula are represented as strings for display.
    gamma_term = "gamma^M"
    # The text defines beta_k, so for k=M we have beta_M
    beta_M_term = "beta_M"
    
    # The full formula for the factor multiplying sqrt(N)
    final_factor_formula = f"{gamma_term} + (1 - {beta_M_term})"
    
    # The definition of beta_M as given in the text
    beta_M_definition = "product_{t=0 to M} (1 - c * delta_t)"
    
    print("The upper-bound for ||B Q_{0, M}||_infinity can be expressed as sqrt(N) times a factor.")
    print("The factor is given by the following formula:")
    print(final_factor_formula)
    
    print("\nWhere the terms are defined as:")
    print(f"  gamma: A contraction coefficient (0 < gamma < 1) related to the ergodicity of the attention matrices P^(t).")
    print(f"  M: The number of layers (time steps) of the graph attentional network.")
    print(f"  {beta_M_term}: Is defined as {beta_M_definition}, where:")
    print(f"    c: A positive constant.")
    print(f"    delta_t: The value ||D^(t) - I_N||_infinity, which measures how much the diagonal matrix D^(t) deviates from the identity matrix at layer t.")
    
    print("\nIn compliance with the instruction to output each 'number' in the final equation, we list the symbolic components of the factor:")
    print(f"Component 1: {gamma_term}")
    print(f"Component 2: 1")
    print(f"Component 3: -{beta_M_term}")

if __name__ == '__main__':
    display_bound_formula()
