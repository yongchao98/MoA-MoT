import math
from fractions import Fraction

def format_state_b(state_dict):
    """
    Formats a dictionary representing a state of ququint B into a readable string.
    Example: { '0': 0.7071, '2': 0.7071 } -> "0.7071 * |0>_B + 0.7071 * |2>_B"
    """
    if not state_dict:
        return "0"
        
    terms = []
    # Sort by basis state for consistent, predictable output format
    for basis, coeff in sorted(state_dict.items()):
        # Ignore terms with coefficients very close to zero
        if abs(coeff) < 1e-9:
            continue
            
        sign = " + "
        if coeff < 0:
            sign = " - "
        
        # For the first term, handle the sign specially
        if not terms:
            if sign == " - ":
                sign = "-"
            else: # sign is " + "
                sign = ""

        # Format the coefficient's absolute value, omitting it if it's 1.0
        val_str = f"{abs(coeff):.4f}"
        if abs(abs(coeff) - 1.0) < 1e-9:
            term_str = f"|{basis}⟩_B"
        else:
            term_str = f"{val_str}*|{basis}⟩_B"

        terms.append(f"{sign}{term_str}")
    
    # If the state is a superposition, enclose in parentheses
    if len(terms) > 1:
        # Join terms and remove leading " + " from the first term
        return "(" + "".join(terms).lstrip(" +") + ")"
    elif terms:
        # For a single basis state, no parentheses needed
        return terms[0].lstrip(" +")
    else:
        return "0"

def solve_ququint_problem():
    """
    Solves the ququint measurement problem as described.
    """
    # 1. Define the gate Q as its action on basis states.
    # The coefficients are 1/sqrt(2).
    c = 1 / math.sqrt(2)
    Q_op = {
        '0': {'1': c, '2': c},
        '1': {'0': c, '3': c},
        '2': {'1': c, '4': c},
        '3': {'0': c, '2': c},
        '4': {'2': c, '3': c}
    }

    # 2. Define the initial entangled state |\Psi>_AB.
    # |\Psi> = 1/sqrt(5) * (|0>_A|0>_B + |1>_A|1>_B + ...)
    # This is represented as a dictionary mapping (A_state, B_state) -> coefficient.
    initial_coeff = 1 / math.sqrt(5)
    psi_ab = { (str(i), str(i)): initial_coeff for i in range(5) }

    # 3. Apply the Q gate to ququint A.
    # The new state is |\Psi'>_AB = (Q tensor I) |\Psi>_AB.
    psi_prime_ab = {}
    for (a_state, b_state), coeff in psi_ab.items():
        # Get the result of Q|a_state>
        q_a_result = Q_op[a_state]
        for new_a_state, q_coeff in q_a_result.items():
            # The new term in the superposition is |new_a_state>|b_state>
            # with coefficient = initial_coeff * q_coeff
            state_tuple = (new_a_state, b_state)
            new_coeff = coeff * q_coeff
            # Add the new coefficient to any existing coefficient for this state
            psi_prime_ab[state_tuple] = psi_prime_ab.get(state_tuple, 0) + new_coeff
            
    # 4. Group the resulting state by the state of ququint A to analyze measurement.
    grouped_by_a = { str(i): {} for i in range(5) }
    for (a_state, b_state), coeff in psi_prime_ab.items():
        grouped_by_a[a_state][b_state] = coeff

    print("After applying gate Q to ququint A, the system is in a new state.")
    print("If we measure ququint A, it will collapse to one of the following five states with a specific probability:\n")

    # 5. Calculate probabilities and post-measurement states for each outcome.
    for i in range(5):
        a_measured_state = str(i)
        # This is the unnormalized state of B when A is measured as |i>
        b_state_unnormalized = grouped_by_a[a_measured_state]
        
        # The probability is the sum of the squared magnitudes of the coefficients.
        prob = sum(c**2 for c in b_state_unnormalized.values())
        prob_fraction = Fraction(prob).limit_denominator()
        
        if prob < 1e-9:
            print(f"Outcome |{a_measured_state}⟩_A: Probability = 0")
            continue
            
        # The normalization constant is the sqrt of the probability.
        norm = math.sqrt(prob)
        
        # The normalized state of B is found by dividing each coefficient by the norm.
        b_state_normalized = {basis: c / norm for basis, c in b_state_unnormalized.items()}
        
        # Format the output for clarity.
        b_state_str = format_state_b(b_state_normalized)

        print(f"Outcome |{a_measured_state}⟩_A:")
        print(f"  Probability: {prob_fraction.numerator}/{prob_fraction.denominator} ({prob:.2f})")
        print(f"  The final system state is: |{a_measured_state}⟩_A ⊗ {b_state_str}\n")

# Run the simulation and print the results
solve_ququint_problem()