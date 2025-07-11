import re

def cofinality_of_alephs(alpha_str):
    """
    Symbolically computes the cofinality of a cardinal kappa = aleph_alpha.
    The input alpha_str is a string representing the ordinal index alpha.
    This function handles basic cases for alpha.
    """
    # Case: alpha = n (a non-negative integer)
    # An ordinal n is a successor (or 0). aleph_n is regular.
    # cf(aleph_n) = aleph_n.
    if alpha_str.isdigit():
        return f"aleph_{alpha_str}"

    # Case: alpha = omega
    # omega is a limit ordinal. aleph_omega is the first singular cardinal.
    # cf(aleph_omega) = aleph_{cf(omega)} = aleph_0.
    if alpha_str == "omega":
        return "aleph_0"

    # Case: alpha = omega + n
    # omega + n is a successor ordinal. aleph_{omega+n} is regular.
    # cf(aleph_{omega+n}) = aleph_{omega+n}.
    match = re.fullmatch(r"omega\+(\d+)", alpha_str)
    if match:
        return f"aleph_{alpha_str}"

    # Case: alpha = omega_k
    # omega_k is a regular cardinal for k>=1. Thus aleph_{omega_k} is regular.
    # cf(aleph_{omega_k}) = aleph_{omega_k}
    match = re.fullmatch(r"omega_(\d+)", alpha_str)
    if match and int(match.group(1)) > 0:
        return f"aleph_{alpha_str}"
        
    # Default case for other unhandled limit ordinals (e.g., omega*2)
    # Returns the symbolic representation of the cofinality.
    return f"aleph_cf({alpha_str})"

# The theoretical answer to the problem is that the largest mu is the cofinality of kappa.
# The code below demonstrates the answer for two representative cases.

print("The largest value for mu is the cofinality of kappa, cf(kappa).")
print("This depends on whether kappa is a regular or singular cardinal.")
print("-" * 30)

# --- Case 1: kappa is a regular cardinal (e.g., aleph_5) ---
kappa1_index = "5"
kappa1_str = f"aleph_{kappa1_index}"
mu1_str = cofinality_of_alephs(kappa1_index)

print(f"Case 1: Let kappa = {kappa1_str} (a regular cardinal).")
# The final equation is mu = aleph_5.
print(f"Then the largest mu is cf({kappa1_str}) = {mu1_str}.")
mu1_index = re.search(r'\d+$', mu1_str).group()
print("The final equation is mu = " + mu1_str)
print("The number in the index of the final equation is:")
print(mu1_index)
print("-" * 30)

# --- Case 2: kappa is a singular cardinal (e.g., aleph_omega) ---
kappa2_index = "omega"
kappa2_str = f"aleph_{kappa2_index}"
mu2_str = cofinality_of_alephs(kappa2_index)

print(f"Case 2: Let kappa = {kappa2_str} (a singular cardinal).")
# The final equation is mu = aleph_0.
print(f"Then the largest mu is cf({kappa2_str}) = {mu2_str}.")
mu2_index = re.search(r'\d+$', mu2_str).group()
print("The final equation is mu = " + mu2_str)
print("The number in the index of the final equation is:")
print(mu2_index)
print("-" * 30)
