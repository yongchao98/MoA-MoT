import math

def solve_crystal_problem():
    """
    This script solves the two parts of the condensed matter physics problem.
    A. It calculates the dimension of the fiber of the bundle pi.
    B. It calculates the number of independent coefficients for the energy functional E
       using group theory and character tables for the cubic group O.
    """

    # Part A: Dimension of the fiber of pi
    # The field C is a tensor C^i_{j,mu} where i, j, mu each run from 1 to 3.
    # The dimension is the total number of components.
    dim_i = 3
    dim_j = 3
    dim_mu = 3
    fiber_dimension = dim_i * dim_j * dim_mu

    # Part B: Number of coefficients for the energy E
    # This requires a group theory calculation for the cubic group O.

    # Character table information for the cubic rotation group O
    # Classes are: E, 8C3, 3C2 (or 3C4^2), 6C4, 6C2' (or 6C2d)
    group_order = 24
    class_info = {
        'E':   {'order': 1, 'T1_char': 3, 'g2_class': 'E'},
        'C3':  {'order': 8, 'T1_char': 0, 'g2_class': 'C3'},
        'C2':  {'order': 3, 'T1_char': -1, 'g2_class': 'E'},
        'C4':  {'order': 6, 'T1_char': 1, 'g2_class': 'C2'},
        'C2d': {'order': 6, 'T1_char': -1, 'g2_class': 'E'}
    }
    class_names = list(class_info.keys()) # ['E', 'C3', 'C2', 'C4', 'C2d']

    # The dislocation density tensor alpha lives in the representation V = T1 x T1 x T1.
    # Calculate the character of V, chi_V(g) = (chi_T1(g))^3.
    chi_V = {name: info['T1_char']**3 for name, info in class_info.items()}

    # Calculate chi_V(g^2) for each class g.
    chi_V_g2 = {name: chi_V[info['g2_class']] for name, info in class_info.items()}

    # Calculate the character of the symmetric square representation, Sym^2(V).
    # chi_sym2_V(g) = 0.5 * (chi_V(g)^2 + chi_V(g^2))
    chi_sym2_V = {name: 0.5 * (chi_V[name]**2 + chi_V_g2[name]) for name in class_names}

    # The number of coefficients is the number of times the A1 irrep appears in Sym^2(V).
    # This is calculated by the inner product <chi_A1, chi_sym2_V>.
    # Since chi_A1 is 1 for all classes, this simplifies to (1/|G|) * sum_g |g|*chi_sym2_V(g).
    dot_product_sum = sum(class_info[name]['order'] * chi_sym2_V[name] for name in class_names)
    num_coeffs = dot_product_sum / group_order

    # The problem asks to output the numbers in the final equation.
    # Here, we construct a string showing the calculation for the number of coefficients.
    print("Step-by-step calculation for the number of coefficients:")
    
    # Generate the string for the sum terms |g| * chi_sym2_V(g)
    terms_str_list = []
    for name in class_names:
        order = class_info[name]['order']
        chi_val = chi_sym2_V[name]
        terms_str_list.append(f"{order}*{int(chi_val)}")
    
    equation_str = f"NumCoeffs = (1/{group_order}) * [ " + " + ".join(terms_str_list) + " ]"
    result_str = f"            = (1/{group_order}) * {int(dot_product_sum)} = {int(num_coeffs)}"
    
    print(equation_str)
    print(result_str)
    print("-" * 20)

    # Print the final answer in the required "A B" format.
    final_answer = f"{fiber_dimension} {int(num_coeffs)}"
    print("Final Answer:")
    print(final_answer)


# Execute the function to solve the problem and print the results.
solve_crystal_problem()