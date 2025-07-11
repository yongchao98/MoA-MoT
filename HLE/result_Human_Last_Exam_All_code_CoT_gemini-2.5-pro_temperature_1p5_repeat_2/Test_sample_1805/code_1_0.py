def display_nabla_q_Tn(n):
    """
    Calculates and displays the formula for the q-difference quotient of T^n.

    The function also breaks down the resulting equation into its constituent numerical parts
    as requested by the user prompt.

    Args:
        n (int): A non-negative integer representing the power of T.
    """
    print(f"Analysis for n = {n}\n")
    
    if not isinstance(n, int) or n < 0:
        print("Error: n must be a non-negative integer.")
        return

    # --- Part 1: Mathematical Derivation ---
    print("--- Derivation of the Formula ---")
    print("The q-difference quotient (q-derivative) is defined as:")
    print("nabla_q(f(T)) = (f(qT) - f(T)) / (qT - T)")
    print("\nApplying this definition to f(T) = T^n:")
    print(f"nabla_q(T^{n}) = ((qT)^{n} - T^{n}) / (qT - T)")
    print(f"             = (q^{n}*T^{n} - T^{n}) / ((q - 1)T)")
    print(f"             = ((q^{n} - 1) / (q - 1)) * T^{n-1}")
    print(f"The term (q^{n} - 1)/(q - 1) is the q-integer, denoted [n]_q.")
    print(f"[n]_q expands to the polynomial: 1 + q + q^2 + ... + q^{n-1}")
    print(f"So, the final general formula is: nabla_q(T^{n}) = [n]_q * T^{n-1}\n")

    # --- Part 2: Final Equation for the specific n ---
    print("--- Final Equation ---")

    # Handle special cases for n=0 and n=1
    if n == 0:
        equation = f"nabla_q(T^0) = 0"
        power_of_T_str = "N/A"
        q_poly_string = "0"
        coefficients = [0]
        powers_of_q = []
    elif n == 1:
        equation = f"nabla_q(T^1) = 1"
        power_of_T_str = "0 (since T^0 = 1)"
        q_poly_string = "1"
        coefficients = [1]
        powers_of_q = [0]
    else: # General case for n > 1
        # Build the string for the q-integer [n]_q
        q_poly_terms = ["1"]
        if n > 1:
            q_poly_terms.append("q")
        if n > 2:
            q_poly_terms.extend([f"q^{i}" for i in range(2, n)])
        
        q_poly_string = " + ".join(q_poly_terms)

        # Build the T term string
        power_of_T = n - 1
        power_of_T_str = str(power_of_T)
        if power_of_T == 1:
            t_term = "T"
        else:
            t_term = f"T^{power_of_T}"
            
        equation = f"nabla_q(T^{n}) = ({q_poly_string}) * {t_term}"
        
        # Data for the "numbers" part
        coefficients = [1] * n
        powers_of_q = list(range(n))

    print(equation)

    # --- Part 3: Outputting each number in the equation ---
    print("\n--- Details of Each Number in the Equation ---")
    print(f"Initial power 'n' of T: {n}")
    print(f"Resulting power of T: {power_of_T_str}")
    print(f"The q-integer polynomial is: [n]_q = {q_poly_string}")
    print(f"  - Number of terms: {len(coefficients)}")
    print(f"  - Coefficients of the polynomial in q: {coefficients}")
    print(f"  - Powers of q in the polynomial: {powers_of_q}")


if __name__ == '__main__':
    # We use a fixed integer n=5 for a clear, runnable example.
    n_value = 5
    display_nabla_q_Tn(n_value)