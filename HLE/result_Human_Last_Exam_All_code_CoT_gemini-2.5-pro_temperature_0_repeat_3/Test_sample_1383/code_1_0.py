from fractions import Fraction

def solve_growth_rate():
    """
    Calculates and prints the achieved growth rate (W) and the decrease in growth rate (ΔW).
    """
    # True probabilities of each bike winning
    p_true = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    # Incorrectly believed probabilities
    q_incorrect = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    # Decimal odds from the bookmakers (k-for-1)
    odds = [4, 3, 7, 7]

    def format_term(coeff, value):
        """Formats a single term of the equation, e.g., (1/4)*ln(2)."""
        if coeff == 0:
            return None
        
        # Format coefficient string
        coeff_str = f"({coeff})" if coeff.denominator != 1 else str(coeff.numerator)
        if coeff == 1:
            coeff_str = ""
        elif coeff == -1:
            coeff_str = "-"
        else:
            coeff_str = f"({coeff})*"

        # Format value string
        value_str = f"({value})" if value.denominator != 1 else str(value.numerator)
        
        return f"{coeff_str}ln{value_str}"

    def build_equation(terms_dict):
        """Builds the final equation string from a dictionary of terms."""
        parts = []
        # Sort by value for consistent output
        for value in sorted(terms_dict.keys()):
            coeff = terms_dict[value]
            if coeff != 0:
                term_str = format_term(coeff, value)
                parts.append(term_str)
        
        if not parts:
            return "0"
            
        # Join parts, cleaning up signs
        return " + ".join(parts).replace("+ -", "- ").replace(" + (", " + ").replace("*ln", "*ln")


    # --- Calculate W (achieved growth rate with incorrect beliefs) ---
    # W = Σ p_i * ln(q_i * d_i)
    w_terms = {}
    print("Calculating W = Σ p_i * ln(q_i * d_i):")
    for i in range(len(p_true)):
        p = p_true[i]
        q = q_incorrect[i]
        d = odds[i]
        value = q * d
        if value == 1: # ln(1) is 0, so term is 0
            print(f"  Term {i+1}: p_{i+1}*ln(q_{i+1}*d_{i+1}) = ({p})*ln(({q})*{d}) = ({p})*ln(1) = 0")
            continue
        
        print(f"  Term {i+1}: p_{i+1}*ln(q_{i+1}*d_{i+1}) = ({p})*ln(({q})*{d}) = ({p})*ln({value})")
        w_terms[value] = w_terms.get(value, 0) + p
    
    w_equation = build_equation(w_terms)
    print("\nThe achieved doubling rate is:")
    print(f"W = {w_equation}")

    # --- Calculate ΔW (decrease in growth rate) ---
    # ΔW = Σ p_i * ln(p_i / q_i)
    delta_w_terms = {}
    print("\nCalculating ΔW = Σ p_i * ln(p_i / q_i):")
    for i in range(len(p_true)):
        p = p_true[i]
        q = q_incorrect[i]
        value = p / q
        if value == 1: # ln(1) is 0, so term is 0
            print(f"  Term {i+1}: p_{i+1}*ln(p_{i+1}/q_{i+1}) = ({p})*ln(({p})/({q})) = ({p})*ln(1) = 0")
            continue
        
        print(f"  Term {i+1}: p_{i+1}*ln(p_{i+1}/q_{i+1}) = ({p})*ln(({p})/({q})) = ({p})*ln({value})")
        
        # Simplify ln(1/x) = -ln(x) to combine terms
        if value < 1:
            inv_value = 1 / value
            delta_w_terms[inv_value] = delta_w_terms.get(inv_value, 0) - p
        else:
            delta_w_terms[value] = delta_w_terms.get(value, 0) + p

    delta_w_equation = build_equation(delta_w_terms)
    print("\nThe decrease in the doubling rate is:")
    print(f"ΔW = {delta_w_equation}")

solve_growth_rate()