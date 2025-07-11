import math

# Titan architecture constraints
MAX_VAL = 15
MAX_TERMS = 10

def represent_as_expression(n):
    """
    Represents an integer n as an expression of terms <= MAX_VAL.
    The representation used is q*15 + r, which is a list of q 15s and one r.
    """
    if n <= MAX_VAL:
        return [n]
    
    num_fifteens = n // MAX_VAL
    remainder = n % MAX_VAL
    
    expression = [MAX_VAL] * num_fifteens
    if remainder > 0:
        expression.append(remainder)
    return expression

def multiply_expressions(expr1, expr2):
    """
    Simulates the multiplication of two expressions on Titan.
    The resulting expression is the sum of the representations of all cross-products.
    """
    result_expr = []
    for term1 in expr1:
        for term2 in expr2:
            product = term1 * term2
            product_expr = represent_as_expression(product)
            # Check for term limit violation at each sub-step
            if len(product_expr) > MAX_TERMS:
                print(f"  - Sub-calculation failed: Multiplying {term1} by {term2} gives {product}.")
                print(f"  - The resulting expression {product_expr} has {len(product_expr)} terms.")
                raise ValueError("Expression exceeds 10-term limit during sub-calculation.")
            result_expr.extend(product_expr)
    
    # Check for final term limit violation
    if len(result_expr) > MAX_TERMS:
        print(f"  - Final expression has {len(result_expr)} terms.")
        raise ValueError("Expression exceeds 10-term limit.")
        
    return result_expr

def main():
    """
    Main function to trace the feasibility of the calculation.
    We will trace the calculation of the numerical coefficient F.
    A crude approximation F = 8 * 7 * 12 is sufficient to show the failure.
    """
    print("Attempting to calculate the escape velocity of Pandora on Titan.")
    print("This requires calculating a numerical coefficient, let's use a simplified one: F = 8 * 7 * 12.")
    print(f"Constraint: 4-bit integers (0-{MAX_VAL}), max 10 terms per expression.\n")

    try:
        # Step 1: MOV AX, 8
        ax = represent_as_expression(8)
        print(f"Step 1: MOV AX, 8")
        print(f"AX holds the expression for 8: {ax} ({len(ax)} term(s))\n")

        # Step 2: MUL AX, 7
        print(f"Step 2: MUL AX, 7 (calculating 8 * 7)")
        val_b = 7
        expr_b = represent_as_expression(val_b)
        ax = multiply_expressions(ax, expr_b)
        print(f"AX now holds the expression for 56: {ax} ({len(ax)} term(s))\n")

        # Step 3: MUL AX, 12
        print(f"Step 3: MUL AX, 12 (calculating 56 * 12)")
        val_c = 12
        expr_c = represent_as_expression(val_c)
        
        # This step will fail.
        # The logic is inside multiply_expressions, which will raise an exception.
        # We are multiplying (15+15+15+11) by 12.
        # The first operation will be 15 * 12.
        print("The operation expands to (15*12) + (15*12) + (15*12) + (11*12).")
        print("Let's analyze the first term of this expansion: 15 * 12.")
        ax = multiply_expressions(ax, expr_c)
        
        # This part will not be reached
        print(f"AX now holds the expression for 672: {ax} ({len(ax)} term(s))\n")
        print("Conclusion: The calculation is possible.")
        
    except ValueError as e:
        print(f"\nERROR: Calculation failed. {e}")
        print("\nConclusion: The calculation is NOT possible because intermediate operations generate expressions that exceed the 10-term limit.")
        print("The Titan architecture cannot handle the magnitude of the numbers involved in this physics problem.")

if __name__ == "__main__":
    main()
    print("\n<<<N0>>>")
