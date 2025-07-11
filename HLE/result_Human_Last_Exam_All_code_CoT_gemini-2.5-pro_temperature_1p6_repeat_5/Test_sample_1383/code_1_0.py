import math
from fractions import Fraction

def format_expression(coeffs, names, var_name):
    """
    Formats an expression like a*ln(b) + c*ln(d) into a human-readable string.
    
    Args:
        coeffs (list): A list of Fraction objects representing the coefficients.
        names (list): A list of strings for the variable names (e.g., "ln(2)").
        var_name (str): The name of the quantity being represented (e.g., "W").
        
    Returns:
        str: The formatted string expression.
    """
    terms = list(zip(coeffs, names))
    pos_terms = sorted([t for t in terms if t[0] > 0], key=lambda item: item[1])
    neg_terms = sorted([t for t in terms if t[0] < 0], key=lambda item: item[1])
    
    parts = []
    
    # Process positive terms first
    for coeff, name in pos_terms:
        if coeff == 1:
            term = name
        else:
            term = f"{coeff} * {name}"
        parts.append(term)
        
    # Process negative terms
    for coeff, name in neg_terms:
        num = abs(coeff)
        if num == 1:
            term = f"- {name}"
        else:
            term = f"- {num} * {name}"
        parts.append(term)

    if not parts:
        return f"{var_name} = 0"
    
    # Join the parts into a final string
    final_str = f"{var_name} = {parts[0]}"
    for part in parts[1:]:
        if part.startswith("-"):
            final_str += f" {part}"
        else:
            final_str += f" + {part}"
            
    return final_str

def solve():
    """
    Solves the betting problem and prints the results.
    """
    # True probabilities, incorrect probabilities, and odds
    p = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
    q = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 8), Fraction(1, 8)]
    o = [4, 3, 7, 7]

    # --- Step 1: Optimal Strategy and Growth Rate W* ---
    # Under true probabilities p, only bike 1 has positive expectation (p1*o1 = 1/2*4 = 2 > 1).
    p1, o1 = p[0], o[0]
    b1_star = Fraction((p1 * o1 - 1), (o1 - 1))  # This is 1/3

    # We can express log returns in terms of coefficients for ln(2) and ln(3).
    # ln(a/b) = ln(a) - ln(b), ln(x^y) = y*ln(x)
    # Return on win: 1 - b1_star + b1_star * o1 = 1 - 1/3 + 4/3 = 2. ln(2) -> (ln(2): 1, ln(3): 0)
    # Return on loss: 1 - b1_star = 1 - 1/3 = 2/3. ln(2/3) -> (ln(2): 1, ln(3): -1)
    
    w_star_ln2 = p1 * 1 + (1 - p1) * 1
    w_star_ln3 = p1 * 0 + (1 - p1) * (-1)

    # --- Step 2: Suboptimal Strategy and Achieved Growth Rate W ---
    # Using incorrect probabilities q, only bike 2 seems favorable (q2*o2 = 1/2*3 = 1.5 > 1).
    q2, o2 = q[1], o[1]
    b2_sub = Fraction((q2 * o2 - 1), (o2 - 1))  # This is 1/4

    # Calculate actual growth rate W using true probabilities p and bet b2_sub.
    p2 = p[1]
    # Return if bike 2 wins (prob p2): 1 - b2_sub + b2_sub*o2 = 1-1/4+3/4=3/2. ln(3/2)->(ln(2):-1,ln(3):1)
    # Return if bike 2 loses (prob 1-p2): 1 - b2_sub = 3/4. ln(3/4)->(ln(2):-2, ln(3):1)
    
    w_ln2 = p2 * (-1) + (1 - p2) * (-2)
    w_ln3 = p2 * 1 + (1 - p2) * 1

    # --- Step 3: Calculate the Difference ΔW ---
    delta_w_ln2 = w_star_ln2 - w_ln2
    delta_w_ln3 = w_star_ln3 - w_ln3

    # --- Step 4: Print the formatted results ---
    # The question asks for the doubling rate W and the decrease ΔW.
    
    # Calculate W expression
    w_expression = format_expression([w_ln2, w_ln3], ["ln(2)", "ln(3)"], "W")
    print(f"The achieved doubling rate is:\n{w_expression}")

    # Calculate ΔW expression
    delta_w_expression = format_expression([delta_w_ln2, delta_w_ln3], ["ln(2)", "ln(3)"], "ΔW")
    print(f"\nThe decrease in doubling rate is:\n{delta_w_expression}")

solve()
<<<W = ln(3) - 7/4 * ln(2); ΔW = 11/4 * ln(2) - 3/2 * ln(3)>>>