def solve_moment_problem():
    """
    This script explains the proof for the following problem:

    Suppose f is a Schwartz class function from R to R such that the integral
    of x^k * f(x) over R is 0 for all non-negative integers k.
    Does it follow that f(x) = 0 for all x?

    The answer is Yes. This script lays out the argument.
    """
    print("--- The Moment Problem for Schwartz Functions ---")
    print("\nProblem: Given a Schwartz function f(x) where all its moments are zero, is f(x) necessarily the zero function?")
    print("The k-th moment is defined as M_k = integral(x^k * f(x) dx) from -inf to +inf.")
    print("Condition: M_k = 0 for all k = 0, 1, 2, ...")
    print("-" * 50)

    print("\nStep 1: The Fourier Transform and its Derivatives")
    print("Let F(xi) denote the Fourier transform of f(x).")
    print("F(xi) = integral(f(x) * exp(-2*pi*i*x*xi) dx)")
    print("\nA key property of the Fourier transform is that the derivatives of F(xi) are related to the moments of f(x).")
    print("By differentiating F(xi) k times with respect to xi, we get:")
    print("F^(k)(xi) = d^k/d(xi)^k [F(xi)] = integral(f(x) * (-2*pi*i*x)^k * exp(-2*pi*i*x*xi) dx)")
    print("(Note: Differentiating under the integral is allowed because f(x) is a Schwartz function.)")
    print("-" * 50)

    print("\nStep 2: Connecting to the Moments")
    print("Let's evaluate the k-th derivative of F(xi) at the point xi = 0.")
    print("F^(k)(0) = integral(f(x) * (-2*pi*i*x)^k * exp(0) dx)")
    print("F^(k)(0) = (-2*pi*i)^k * integral(x^k * f(x) dx)")
    print("The term 'integral(x^k * f(x) dx)' is precisely the k-th moment of f(x), M_k.")
    print("So, we have the relation: F^(k)(0) = (-2*pi*i)^k * M_k")
    print("-" * 50)

    print("\nStep 3: Applying the Given Condition")
    print("The problem states that all moments M_k are zero.")
    print("M_k = integral(x^k * f(x) dx) = 0 for all k = 0, 1, 2, ...")
    print("\nSubstituting this into our equation for the derivatives of F at the origin:")
    moment = 0
    equation = f"F^(k)(0) = (-2*pi*i)^k * {moment}"
    print(equation)
    final_derivative_value = 0
    print(f"F^(k)(0) = {final_derivative_value} for all k = 0, 1, 2, ...")
    print("\nThis means that the Fourier transform F(xi) and all of its derivatives are zero at xi=0.")
    print("-" * 50)

    print("\nStep 4: Using Properties of Schwartz Functions")
    print("A crucial theorem in Fourier analysis states that the Fourier transform of a Schwartz function is not only infinitely differentiable but also analytic.")
    print("An analytic function is one that is locally equal to its Taylor series.")
    print("\nThe Taylor series of F(xi) expanded around xi = 0 is:")
    print("Series = F(0) + F'(0)*xi + (F''(0)/2!)*xi^2 + ... = sum_k (F^(k)(0)/k!) * xi^k")
    print("\nSince we have shown that F^(k)(0) = 0 for all k, every coefficient in the Taylor series is zero.")
    print("Series = 0 + 0*xi + 0*xi^2 + ... = 0")
    print("\nBecause F(xi) is analytic, it must be equal to its Taylor series. Therefore, F(xi) must be identically zero for all xi.")
    print("-" * 50)

    print("\nStep 5: Conclusion via Fourier Inversion")
    print("The Fourier Inversion Theorem allows us to recover f(x) from its transform F(xi):")
    print("f(x) = integral(F(xi) * exp(2*pi*i*x*xi) d(xi))")
    print("\nSince we have concluded that F(xi) = 0 for all xi, the integral becomes trivial:")
    final_integral = "integral(0 * exp(2*pi*i*x*xi) d(xi))"
    final_value = 0
    print(f"f(x) = {final_integral} = {final_value}")
    print("\nThis shows that f(x) must be the zero function.")
    print("-" * 50)

    print("\nFinal Answer:")
    print("Yes, if a Schwartz function f(x) has all its moments equal to zero, it follows that f(x) = 0.")
    print("The final equation is:")
    print("f(x) = 0")

if __name__ == '__main__':
    solve_moment_problem()
