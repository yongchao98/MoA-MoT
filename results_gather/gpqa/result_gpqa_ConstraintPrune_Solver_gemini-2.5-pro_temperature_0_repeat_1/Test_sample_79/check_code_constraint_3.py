import math
from sympy.solvers.diophantine.diophantine import diop_linear
from sympy import symbols, sympify

def is_prime(n):
    """Checks if a number is prime."""
    # Ensure n is a standard integer for primality test
    try:
        n_int = int(n)
    except (ValueError, TypeError):
        return False
        
    if n_int != n or n_int <= 1:
        return False
    if n_int <= 3:
        return True
    if n_int % 2 == 0 or n_int % 3 == 0:
        return False
    i = 5
    while i * i <= n_int:
        if n_int % i == 0 or n_int % (i + 2) == 0:
            return False
        i += 6
    return True

def check_correctness():
    """
    Checks the correctness of the LLM's answer by systematically finding all
    possible integer solutions for each option and evaluating them against the
    "maximize primes" hypothesis.
    """
    # The target string is ACAGTGACC
    # Counts: A=3, C=3, G=2, T=1
    # Target value V = 3*v_A + 3*v_C + 2*v_G + v_T

    # From the examples, we have two constraints:
    # 1) v_A = 115 - 2*v_G
    # 2) v_C = 61 - 2*v_T

    # Substituting these into the expression for V gives a Diophantine equation for each option V:
    # 4*v_G + 5*v_T = 528 - V

    options = {'A': 351, 'B': 333, 'C': 315, 'D': 185}
    results = {}
    
    v_G_sym, v_T_sym, t_sym = symbols('v_G v_T t', integer=True)

    for option_name, V in options.items():
        c = 528 - V
        
        # Solve the Diophantine equation 4*v_G + 5*v_T = c
        # The equation passed to diop_linear is ax+by+c=0 form
        sol = diop_linear(4*v_G_sym + 5*v_T_sym - c)
        if not sol:
            results[option_name] = 0
            continue
        
        v_G_expr, v_T_expr = sol[0] if isinstance(sol, set) else sol

        max_primes_for_option = 0
        
        # Iterate through a reasonable range for the parameter 't' to find all
        # solutions where v_A, v_C, v_G, v_T are positive integers.
        # A range of -200 to 200 is more than sufficient.
        for t in range(-200, 201):
            v_G = v_G_expr.subs(t_sym, t)
            v_T = v_T_expr.subs(t_sym, t)
            
            v_A = 115 - 2*v_G
            v_C = 61 - 2*v_T
            
            values = [v_A, v_C, v_G, v_T]
            if all(isinstance(v, sympify(1).__class__) and v > 0 for v in values):
                prime_count = sum(1 for v in values if is_prime(v))
                if prime_count > max_primes_for_option:
                    max_primes_for_option = prime_count
        
        results[option_name] = max_primes_for_option

    # The LLM's reasoning is that option C is the *unique* answer that maximizes the prime count.
    # Let's check this claim.
    
    llm_answer_option = 'C'
    
    # Find the maximum prime count across all options
    if not results or max(results.values()) == 0:
         return "Could not find any valid positive integer solutions for any option."
         
    max_overall_primes = max(results.values())
    
    # Find all options that achieve this maximum
    options_with_max_primes = [opt for opt, count in results.items() if count == max_overall_primes]
    
    if results[llm_answer_option] < max_overall_primes:
        return (f"Incorrect. The answer selects option {llm_answer_option}, which yields a maximum of "
                f"{results[llm_answer_option]} primes. However, option(s) "
                f"{', '.join(options_with_max_primes)} yield(s) a higher count of {max_overall_primes} primes.")

    if len(options_with_max_primes) > 1:
        return (f"Incorrect. The reasoning provided in the answer is flawed. It claims that option C *uniquely* "
                f"maximizes the number of primes with a count of 3. However, a full analysis shows that the "
                f"maximum number of primes is indeed {max_overall_primes}, but this is achieved by multiple options: "
                f"{', '.join(options_with_max_primes)}. For example, for option A (351), the values {{v_A=89, v_C=11, v_G=13, v_T=25}} "
                f"also yield 3 primes. Since the 'maximize primes' constraint does not lead to a unique choice, "
                f"the reasoning is invalid.")

    # This case means the LLM's chosen option is the unique one with the max primes.
    return "Correct"

# Run the check and print the result.
result = check_correctness()
print(result)