import sympy

def check_correctness():
    """
    This function checks the correctness of the provided answer by symbolically
    deriving the relationship using the sympy library.

    The problem involves transforming a number density function from parallax (plx)
    space to distance (r) space.

    1.  **Given:** The number of stars per unit parallax, dN/d(plx), is
        proportional to 1/plx^5. This is the standard physical interpretation
        of such a statement.
    2.  **Relationship:** Parallax is inversely proportional to distance, plx ∝ 1/r.
    3.  **Goal:** Find the number of stars per unit distance, dN/dr.
    4.  **Method:** Use the chain rule for transforming distributions:
        dN/dr = (dN/d(plx)) * |d(plx)/dr|
    """
    try:
        # Define symbols for distance (r) and constants of proportionality (k, C).
        # We assume distance and constants are positive, as is physically necessary.
        r = sympy.Symbol('r', positive=True)
        k = sympy.Symbol('k', positive=True) # For plx = k/r
        C = sympy.Symbol('C', positive=True) # For dN/d(plx) = C/plx^5

        # --- Step 1: Define the fundamental relationships ---
        # Parallax (plx) is inversely proportional to distance (r)
        plx = k / r

        # The number of stars per unit parallax, n_plx, is proportional to 1/plx^5
        n_plx = C / plx**5

        # --- Step 2: Calculate the Jacobian of the transformation ---
        # We need the absolute value of the derivative of plx with respect to r.
        d_plx_dr = sympy.diff(plx, r)
        abs_jacobian = sympy.Abs(d_plx_dr)

        # --- Step 3: Apply the transformation rule ---
        # The number of stars per unit distance, n_r, is given by:
        # n_r = n_plx(r) * |d(plx)/dr|
        # First, express n_plx in terms of r by substituting plx = k/r
        n_plx_in_terms_of_r = n_plx.subs(plx, k / r)

        # Now, calculate n_r by multiplying by the Jacobian
        n_r = n_plx_in_terms_of_r * abs_jacobian

        # --- Step 4: Simplify the final expression ---
        simplified_n_r = sympy.simplify(n_r)

        # --- Step 5: Check the proportionality ---
        # The final answer claims the proportionality is r^3.
        # To check this, we divide our result by r^3. If the result is a constant
        # (i.e., not dependent on r), then the proportionality is correct.
        proportionality_check = simplified_n_r / (r**3)

        # is_constant() checks if the expression is free of the variable 'r'
        if proportionality_check.is_constant():
            # The derivation confirms the result is proportional to r^3.
            # The provided answer is 'C', which corresponds to ~r^3.
            # The logic in the final answer is sound.
            return "Correct"
        else:
            # If the proportionality is not r^3, the answer is incorrect.
            # We can determine the actual power of r from our derivation.
            power_of_r = sympy.degree(simplified_n_r.as_poly(r))
            return (f"Incorrect. The symbolic derivation shows that the number of stars "
                    f"per unit distance (dN/dr) is proportional to r^{power_of_r}, "
                    f"not r^3. The provided answer 'C' is wrong.")

    except Exception as e:
        return f"An error occurred during the symbolic check: {e}"

# Run the check and print the result
result = check_correctness()
print(result)