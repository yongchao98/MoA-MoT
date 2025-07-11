def compute_integral_lambda_321_on_M3():
    """
    This function calculates the integral of the product of lambda classes
    lambda_3 * lambda_2 * lambda_1 on the moduli space of stable curves
    of genus 3, denoted M_3_bar.
    """

    # The genus of the curves is g=3.
    g = 3

    # The dimension of the moduli space M_g_bar is 3*g - 3.
    dim_space = 3 * g - 3

    # The integrand is the product of lambda classes lambda_c1, lambda_c2, lambda_c3.
    # The class lambda_i has degree i.
    class_indices = [3, 2, 1]
    integrand_degree = sum(class_indices)

    print(f"The integral to compute is ∫_{{M̄_{g}}} λ_{class_indices[0]}λ_{class_indices[1]}λ_{class_indices[2]} for g={g}.")
    print(f"The dimension of the space M̄_{g} is {dim_space}.")
    print(f"The degree of the integrand is {class_indices[0]} + {class_indices[1]} + {class_indices[2]} = {integrand_degree}.")
    print("The degree matches the dimension, so the integral is a number.")

    # We use a fundamental result from intersection theory on moduli spaces.
    # A theorem by C. Faber and R. Pandharipande (2000) states that
    # for any genus g >= 3, the product of classes λ_g * λ_{g-1} is zero
    # in the Chow ring of the compactified moduli space M_g_bar.
    # For g = 3, this means λ_3 * λ_2 = 0.
    
    print("\nApplying the Faber-Pandharipande theorem for g=3: λ_3 * λ_2 = 0.")

    # We can now simplify the integrand:
    # λ_3 * λ_2 * λ_1 = (λ_3 * λ_2) * λ_1 = 0 * λ_1 = 0
    
    print("The integrand simplifies as follows, showing each number in the equation:")
    print(f"λ_{class_indices[0]} * λ_{class_indices[1]} * λ_{class_indices[2]} = (λ_{class_indices[0]} * λ_{class_indices[1]}) * λ_{class_indices[2]} = 0 * λ_{class_indices[2]} = 0")
    
    # The integral of the zero class over any space is zero.
    integral_result_num = 0
    integral_result_den = 1

    print("\nThe integral of the zero class is 0.")
    print("\nFinal Result:")
    print(f"{integral_result_num}/{integral_result_den}")

compute_integral_lambda_321_on_M3()