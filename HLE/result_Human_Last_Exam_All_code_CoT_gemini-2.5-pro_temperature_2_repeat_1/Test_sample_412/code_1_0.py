import math

def run_analysis():
    """
    Analyzes the non-standard analysis definitions of continuity.
    """
    print("Analyzing the properties of a continuous function using a non-standard analysis model.\n")

    # --- Conceptual Model ---
    # We'll represent a hyperreal number x as a tuple (st, inf), where
    # 'st' is the standard part (a regular float) and
    # 'inf' is the infinitesimal part (a coefficient for an infinitesimal epsilon ε).
    # A standard number has an infinitesimal part of 0, e.g., 5 is represented as (5.0, 0.0).
    # A non-standard number like 5 + 2ε would be (5.0, 2.0).

    def are_close(x, y):
        """Checks if two hyperreals x and y are infinitesimally close (x ~ y)."""
        # Two numbers are infinitesimally close if their standard parts are equal.
        return x[0] == y[0]

    def print_point(name, p):
        """Helper to print a hyperreal point."""
        if p[1] == 0:
            print(f"{name} = {p[0]} (a standard point)")
        else:
            print(f"{name} = {p[0]} + {p[1]}ε (a non-standard point)")

    # --- Functions to Test ---
    # Continuous function: f(x) = x^2
    def f_continuous(p):
        st, inf = p
        # (st + inf*ε)^2 = st^2 + 2*st*inf*ε + inf^2*ε^2
        # We ignore higher-order infinitesimals like ε^2.
        return (st**2, 2 * st * inf)

    # Discontinuous function: f(x) = floor(x)
    def f_discontinuous(p):
        st, inf = p
        # For a point x = st + inf*ε, if inf is positive, x > st, so floor(x) = st.
        # If inf is negative, x < st, so floor(x) = st - 1.
        if inf >= 0:
            return (math.floor(st), 0.0)
        else:
            return (math.floor(st) - 1, 0.0)

    # --- Analyzing Property B ---
    # Property B: ∀ x₀ ∈ X, ∀ x₁ ∈ X*: x₀ ~ x₁ ⟹ f(x₀) ~ f(x₁)
    print("--- Testing Property B: x₀ ~ x₁ ⟹ f(x₀) ~ f(x₁) ---\n")
    print("This property states that if a non-standard point x₁ is infinitely close to a standard point x₀,")
    print("then their images under f must also be infinitely close.\n")

    # 1. Test with the continuous function f(x) = x^2
    print("1. Testing with a continuous function: f(x) = x² at x₀ = 3")
    x0_std = (3.0, 0.0)       # A standard point x₀
    x1_nonstd = (3.0, 0.1)  # A non-standard point x₁ infinitesimally close to x₀
    
    premise = are_close(x0_std, x1_nonstd)
    
    fx0 = f_continuous(x0_std)
    fx1 = f_continuous(x1_nonstd)
    conclusion = are_close(fx0, fx1)
    
    print_point("x₀", x0_std)
    print_point("x₁", x1_nonstd)
    print(f"Are x₀ and x₁ infinitely close? (x₀ ~ x₁): {premise}")
    
    print_point("f(x₀)", fx0)
    print_point("f(x₁)", fx1)
    print(f"Are f(x₀) and f(x₁) infinitely close? (f(x₀) ~ f(x₁)): {conclusion}")
    
    print(f"Implication (x₀ ~ x₁ ⟹ f(x₀) ~ f(x₁)) holds: {not premise or conclusion}")
    print("Conclusion: Property B holds for the continuous function f(x)=x².\n")

    # 2. Test with the discontinuous function f(x) = floor(x)
    print("2. Testing with a discontinuous function: f(x) = floor(x) at x₀ = 3")
    # A standard point, which is a point of discontinuity for floor()
    x0_std_floor = (3.0, 0.0)
    
    # We choose a non-standard point x₁ just below x₀
    x1_nonstd_floor = (3.0, -0.1) 

    premise_floor = are_close(x0_std_floor, x1_nonstd_floor)
    
    fx0_floor = f_discontinuous(x0_std_floor)
    fx1_floor = f_discontinuous(x1_nonstd_floor)
    conclusion_floor = are_close(fx0_floor, fx1_floor)
    
    print_point("x₀", x0_std_floor)
    print_point("x₁", x1_nonstd_floor)
    print(f"Are x₀ and x₁ infinitely close? (x₀ ~ x₁): {premise_floor}")

    print_point("f(x₀)", fx0_floor)
    print_point("f(x₁)", fx1_floor)
    print(f"Are f(x₀) and f(x₁) infinitely close? (f(x₀) ~ f(x₁)): {conclusion_floor}")

    print(f"Implication (x₀ ~ x₁ ⟹ f(x₀) ~ f(x₁)) holds: {not premise_floor or conclusion_floor}")
    print("Conclusion: Property B fails for the discontinuous function f(x)=floor(x) because a point")
    print("infinitesimally close to 3 produced an image that is not infinitely close to f(3).\n")

    print("--- Final Conclusion ---")
    print("Property B is the only property that correctly distinguishes continuous from discontinuous functions.")
    print("A is trivial, C is too strong (fails for many continuous functions), and D, E, F describe other properties like injectivity or being a closed map, not continuity.")
    print("\nTherefore, B is the correct equivalent property for a map being continuous.")


if __name__ == "__main__":
    run_analysis()