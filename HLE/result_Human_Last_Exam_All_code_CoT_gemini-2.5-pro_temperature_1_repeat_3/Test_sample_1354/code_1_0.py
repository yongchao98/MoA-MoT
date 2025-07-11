import numpy as np

def solve_trace_covariance():
    """
    Calculates the trace of the covariance matrix for the specified sampling procedure.
    
    The problem simplifies to calculating 1 - (E[(a-b)/(a+b)])^2.
    This is derived as follows:
    1. Tr(Cov(v)) = Tr(Cov(d)) because v is an orthogonal transformation of d.
    2. Tr(Cov(d)) = E[||d||^2] - ||E[d]||^2.
    3. The vector d is constructed to have norm ||d||=1, so E[||d||^2]=1.
    4. E[d] simplifies because E[c/||c||]=0 due to the symmetric distribution of c.
       This leaves ||E[d]||^2 = (E[(a-b)/(a+b)])^2.
    5. E[(a-b)/(a+b)] is calculated using the properties of the Beta distribution,
       where a/(a+b) ~ Beta(alpha, beta).
    """

    # --- Given parameters ---
    alpha = 3.0
    beta = 2.0
    
    # 1. Calculate the mean of the Beta-distributed variable U = a/(a+b)
    # For U ~ Beta(alpha, beta), E[U] = alpha / (alpha + beta)
    mean_U = alpha / (alpha + beta)
    
    # 2. Calculate E[(a-b)/(a+b)] = E[2*U - 1] = 2*E[U] - 1
    exp_val_d1 = 2 * mean_U - 1
    
    # 3. Calculate the squared norm of E[d]
    norm_sq_Ed = exp_val_d1**2
    
    # 4. Calculate the final trace: Tr(Cov(v)) = 1 - ||E[d]||^2
    trace_covariance = 1 - norm_sq_Ed
    
    # --- Print the results step-by-step ---
    print(f"The trace of the covariance matrix is analytically solved by the equation: 1 - (E[(a-b)/(a+b)])^2")
    print(f"Step 1: Calculate E[a/(a+b)], which is the mean of a Beta({alpha}, {beta}) distribution.")
    print(f"E[a/(a+b)] = {alpha} / ({alpha} + {beta}) = {mean_U:.2f}")
    
    print("\nStep 2: Calculate E[(a-b)/(a+b)].")
    print(f"E[(a-b)/(a+b)] = 2 * E[a/(a+b)] - 1 = 2 * {mean_U:.2f} - 1 = {exp_val_d1:.2f}")

    print("\nStep 3: Calculate the trace.")
    print(f"Trace = 1 - (E[(a-b)/(a+b)])^2")
    print(f"Trace = 1 - ({exp_val_d1:.2f})^2 = 1 - {norm_sq_Ed:.2f} = {trace_covariance:.2f}")

solve_trace_covariance()
<<<0.96>>>