import sympy as sp

def analyze_unboundedness(s_val, p_val):
    """
    Analyzes the condition for the functional J_t to be unbounded from below by
    comparing the scaling exponents of its terms.

    Args:
        s_val (float): The parameter 's' from the scaling transformation.
        p_val (float): The parameter 'p' from the L^p norm.
    """
    s, p, t = sp.symbols('s p t', real=True, positive=True)

    print(f"--- Analysis for s = {s_val} and p = {p_val} ---")

    # Define exponents symbolically
    kinetic_exponent_x = 2 * s
    potential_exponent = p * (1 + s) / 2 - (s + 1)

    # To show the functional is unbounded below, we can choose a trial function u(x,y)
    # that varies predominantly in x. For such a function, the kinetic energy of u_t
    # will be dominated by the term scaling with t^(2s).
    
    # We compare the potential exponent with this dominant kinetic exponent.
    # J_t -> -inf if potential_exponent > kinetic_exponent_x
    
    ke_exp_val = kinetic_exponent_x.subs(s, s_val)
    pe_exp_val = potential_exponent.subs({s: s_val, p: p_val})
    
    print("For a trial function varying mostly in x, the dominant kinetic term scales as t^exp, where:")
    print(f"Kinetic Exponent = 2*s = 2*{s_val} = {float(ke_exp_val)}")
    
    print("\nThe potential term involving L^p norm scales as t^exp, where:")
    final_eq_str = f"Potential Exponent = p(1+s)/2 - (s+1) = {p_val}(1+{s_val})/2 - ({s_val}+1) = {float(pe_exp_val)}"
    print(final_eq_str)

    is_unbounded = pe_exp_val > ke_exp_val

    print("\nFor J_t to be unbounded below, we need: Potential Exponent > Kinetic Exponent")
    print(f"Checking the inequality: {float(pe_exp_val)} > {float(ke_exp_val)} is {is_unbounded}")

    # The condition in the question is p > 2(1+3s)/(1+s).
    # This is mathematically equivalent to Potential Exponent > Kinetic Exponent.
    
    threshold = 2 * (1 + 3 * s_val) / (1 + s_val)
    print(f"\nThe condition from question (a) is p > 2(1+3s)/(1+s).")
    print(f"For s={s_val}, this means p > {threshold:.4f}.")
    print(f"Our given p={p_val} satisfies this condition: {p_val > threshold}.")
    print("\nThis confirms that if the condition from the question holds, J_t is indeed unbounded below.")

# Example case: s=2, and p=5 which satisfies p > 2(1+3*2)/(1+2) = 4.66...
print("Demonstrating the logic for question (a) with an example.\n")
analyze_unboundedness(s_val=2.0, p_val=5.0)
