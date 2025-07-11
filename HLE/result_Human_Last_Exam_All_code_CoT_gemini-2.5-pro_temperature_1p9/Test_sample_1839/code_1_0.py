import re

def get_semidistributive_mu(kappa_str: str):
    """
    For a forcing notion P with density kappa, this function calculates the 
    largest cardinal mu such that P is necessarily (mu, kappa^+)-semidistributive.

    The answer from set theory is that mu is the cofinality of kappa, denoted cf(kappa).

    This script computes mu for a given kappa on the Aleph scale (e.g., "aleph_5" 
    or "aleph_omega").

    Args:
        kappa_str: A string representing the cardinal kappa.
    """
    
    print(f"The problem asks for the largest cardinal 'mu' such that any forcing notion P with density kappa = {kappa_str} is necessarily (mu, kappa^+)-semidistributive.")
    print("The solution from advanced set theory is that mu = cf(kappa), the cofinality of kappa.")
    print("-" * 20)

    # Regular expression to parse "aleph_index"
    match = re.match(r"aleph_(\w+)", kappa_str)
    if not match:
        print(f"Error: Invalid input format for kappa: '{kappa_str}'.")
        print("Please use the format 'aleph_n' (for integers n) or 'aleph_omega'.")
        return

    index_str = match.group(1)
    
    # In the Aleph scale, for a cardinal kappa = aleph_alpha, its cofinality depends
    # on the cofinality of the ordinal index alpha.
    # We handle two representative cases:
    # 1. kappa is a regular cardinal, e.g., aleph_n for finite n.
    # 2. kappa is a singular cardinal, e.g., aleph_omega.

    mu_str = ""
    kappa_index_repr = ""
    mu_index_repr = ""

    if index_str == "omega":
        # Case: kappa = aleph_omega. This is a singular cardinal.
        # Its cofinality is cf(aleph_omega) = aleph_0.
        kappa_index_repr = "ω"
        mu_index_repr = "0"
        mu_str = f"aleph_{mu_index_repr}"
        
    else:
        try:
            index_val = int(index_str)
            # Case: kappa = aleph_n for a finite number n. These cardinals are regular.
            # Their cofinality is the cardinal itself. cf(aleph_n) = aleph_n.
            kappa_index_repr = str(index_val)
            mu_index_repr = str(index_val)
            mu_str = kappa_str
        except ValueError:
            print(f"Error: Unsupported index '{index_str}'. This script only supports integer indices and 'omega'.")
            return

    # Print the final result in a structured equation format.
    # The format 'aleph_o' is a stand-in for aleph subscript omega.
    final_equation = f"μ = cf(ℵ_{kappa_index_repr}) = ℵ_{mu_index_repr}"
    
    print("Calculation:")
    print(f"For kappa = {kappa_str} (ℵ_{kappa_index_repr}):")
    print(f"μ = cf({kappa_str}) = {mu_str}")

    print("\nFinal Equation:")
    print(final_equation)

# --- Execution of the code ---
# We run an example case for a singular cardinal, which is more illustrative.
# The user can change "aleph_omega" to another value like "aleph_5" to test other cases.
get_semidistributive_mu("aleph_omega")
