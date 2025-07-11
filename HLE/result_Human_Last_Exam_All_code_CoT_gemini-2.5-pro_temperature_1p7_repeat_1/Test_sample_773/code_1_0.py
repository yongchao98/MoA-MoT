import math

def calculate_mass(n, q, q_v):
    """
    Calculates the total mass based on the derived formula.

    The formula for the total mass Q is:
    Q = (q_v / (q_v - 1)) * Product_{j=1 to n-1} (q^j / (q^j - 1))
    """

    if q <= 1 or q_v <= 1 or n < 1:
        print("Parameters must be n >= 1, q > 1, q_v > 1")
        return

    # Calculate the product term: Product_{j=1 to n-1} (q^j / (q^j - 1))
    product_val = 1
    # Term breakdown for printing
    product_terms_str = []
    
    for j in range(1, n):
        term_numerator = q**j
        term_denominator = q**j - 1
        term = term_numerator / term_denominator
        product_val *= term
        product_terms_str.append(f"({term_numerator} / {term_denominator})")

    # Calculate the q_v prefactor
    prefactor_num = q_v
    prefactor_den = q_v - 1
    prefactor = prefactor_num / prefactor_den
    
    # Final result
    total_mass = prefactor * product_val

    # Printing the final equation with all numbers
    print("The final equation for the total mass is:")
    
    equation_str = f"Mass = ({prefactor_num} / {prefactor_den})"
    if product_terms_str:
        equation_str += " * " + " * ".join(product_terms_str)
    
    # Calculate the intermediate results for printing
    intermediate_results = [f"{prefactor}"]
    if product_terms_str:
      intermediate_results.append(f"{product_val}")
      
    full_equation_str = f"{equation_str} = {' * '.join(intermediate_results)} = {total_mass}"

    print(full_equation_str)

    return total_mass

# Since the parameters n, q, and q_v are not provided in the problem,
# we will use the simplest non-trivial case to provide a concrete numerical answer.
# Let n = 2, q = 2.
# A natural choice is to consider the completion at a place where the residue field
# has the same order as the base field of the function ring, so q_v = q.
# Let q_v = 2.

n_val = 2
q_val = 2
q_v_val = 2

final_answer = calculate_mass(n_val, q_val, q_v_val)

# The final answer in the required format
# To avoid printing None if calculate_mass returns nothing on error
if final_answer is not None:
    # Use math.isclose for float comparison to handle potential precision issues
    if math.isclose(final_answer, round(final_answer)):
        print(f"\n<<< {int(final_answer)} >>>")
    else:
        print(f"\n<<< {final_answer} >>>")
