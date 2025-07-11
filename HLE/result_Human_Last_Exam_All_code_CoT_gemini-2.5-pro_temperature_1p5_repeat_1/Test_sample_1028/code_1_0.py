import textwrap

def explain_identifiability_strategies():
    """
    Explains which strategies help mitigate the identifiability issue in birth-death models
    and identifies the one that does not.
    """

    print("### The Unidentifiability Problem in Birth-Death Models\n")
    explanation = (
        "When fitting a birth-death model with time-varying speciation (λ(t)) and "
        "extinction (μ(t)) rates to a phylogeny of only extant species, a fundamental "
        "identifiability issue arises. This means that multiple different combinations of "
        "λ(t) and μ(t) functions can result in the exact same likelihood for the "
        "observed phylogeny. Therefore, we cannot uniquely determine both functions from "
        "the data without additional information or constraints."
    )
    print(textwrap.fill(explanation, width=80))
    print("\n" + "="*80 + "\n")

    print("### Analysis of Proposed Strategies\n")

    # --- Strategies that HELP ---
    print("Strategies that HELP mitigate the issue:\n")

    print("B. Incorporating Priors:")
    explanation_b = (
        "This helps. Informative priors in a Bayesian framework add external information, "
        "constraining the parameters to biologically plausible ranges and making them "
        "identifiable from the posterior distribution."
    )
    print(textwrap.fill(explanation_b, width=80, initial_indent="  ", subsequent_indent="  "))
    print("-" * 20)

    print("D & F. Incorporating Fossils:")
    explanation_df = (
        "This helps. Fossils provide direct data on extinct lineages, breaking the ambiguity "
        "between speciation and extinction that exists when only observing survivors."
    )
    print(textwrap.fill(explanation_df, width=80, initial_indent="  ", subsequent_indent="  "))
    print("-" * 20)
    
    print("E & G. Reparameterization:")
    explanation_eg = (
        "This helps. While λ(t) and μ(t) are not identifiable, combinations like the 'pulled "
        "diversification rate' are. By reformulating the model to estimate these identifiable "
        "parameters, we work around the problem."
    )
    print(textwrap.fill(explanation_eg, width=80, initial_indent="  ", subsequent_indent="  "))
    print("-" * 20)

    # --- Strategy that might weakly help or be insufficient ---
    print("A. Piecewise Constant Model:")
    explanation_a = (
        "This is a form of regularization that simplifies the model. While it is an attempt "
        "to mitigate the issue by reducing complexity, it doesn't solve the core problem "
        "within each time slice and is often insufficient alone."
    )
    print(textwrap.fill(explanation_a, width=80, initial_indent="  ", subsequent_indent="  "))
    
    print("\n" + "="*80 + "\n")

    # --- The strategy that does NOT HELP ---
    print("The strategy that does NOT HELP:\n")
    
    print("C. Piecewise High-Degree Polynomial Model:")
    explanation_c = (
        "This does NOT help. In fact, it makes the problem worse. Using high-degree polynomials "
        "dramatically increases model complexity and flexibility. This provides far more ways "
        "for λ(t) and μ(t) to covary and produce the same likelihood, exacerbating the "
        "identifiability problem. It is the opposite of the required strategies of adding "
        "information or simplifying the model."
    )
    print(textwrap.fill(explanation_c, width=80, initial_indent="  ", subsequent_indent="  "))
    

if __name__ == '__main__':
    explain_identifiability_strategies()
    print("\n<<<C>>>")
