import numpy as np
import matplotlib.pyplot as plt

def smooth_bump_base(t):
    """A smooth function that is positive for t>0 and zero otherwise."""
    # Using a condition to avoid math domain error for np.exp
    # and to ensure the output is zero for t <= 0
    result = np.zeros_like(t, dtype=float)
    mask = t > 0
    result[mask] = np.exp(-1.0 / t[mask])
    return result

def smooth_step(t):
    """A smooth function that transitions from 0 to 1 on the interval [0, 1]."""
    t = np.clip(t, 0, 1) # We only care about the transition in [0,1]
    denom = smooth_bump_base(t) + smooth_bump_base(1 - t)
    # Avoid division by zero at t=0 and t=1, where denom is 0.
    # At t=0, num=0. At t=1, num = denom, so ratio is 1.
    result = np.zeros_like(t, dtype=float)
    mask = (t > 0) & (t < 1)
    result[mask] = smooth_bump_base(t[mask]) / denom[mask]
    result[t >= 1] = 1.0
    return result

def construct_x(t):
    """
    Constructs a smooth function x(t) that covers all real numbers
    and is zero on [-1, 1], enabling y(t)=|x(t)| to be smooth.
    """
    # Create the curve components, shifting and scaling the smooth_step function
    # x(t) will trace the positive reals for t > 1
    # x(t) will trace the negative reals for t < -1
    # x(t) will be 0 for t in [-1, 1]
    
    # Positive part for t > 1
    # We want x(t) to behave like t-2 for large t. Let's make a transition from 1 to 2.
    psi_pos = smooth_step(t - 1)
    x_pos = (t - 2) * psi_pos + psi_pos # Smoothly increases from 0
    
    # Negative part for t < -1
    # We want x(t) to behave like t+2 for large negative t. Transition from -2 to -1.
    psi_neg = smooth_step(-t - 1)
    x_neg = (t + 2) * psi_neg - psi_neg # Smoothly decreases from 0
    
    return x_pos + x_neg

def main():
    """
    Main function to generate and plot the smooth curve whose image is y=|x|.
    This demonstrates that statement B is true.
    """
    # Generate time values
    t = np.linspace(-4, 4, 1000)
    
    # Calculate the components of the curve gamma(t) = (x(t), |x(t)|)
    x_vals = construct_x(t)
    y_vals = np.abs(x_vals)

    # Plotting
    plt.figure(figsize=(8, 8))
    plt.plot(x_vals, y_vals, label=r'$\gamma(t) = (x(t), |x(t)|)$')
    # For reference, plot y=|x| directly
    x_ref = np.linspace(np.min(x_vals), np.max(x_vals), 100)
    y_ref = np.abs(x_ref)
    plt.plot(x_ref, y_ref, 'r--', label=r'$y = |x|$', alpha=0.7)
    
    plt.title('A Smooth Curve whose Image is the Graph of y = |x|')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.axhline(0, color='black',linewidth=0.5)
    plt.axvline(0, color='black',linewidth=0.5)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.legend()
    
    print("This plot shows a smooth curve (in blue) whose image perfectly traces the V-shape of y=|x| (dashed red).")
    print("This visualization confirms that Statement B is True.")
    
    plt.show()

if __name__ == '__main__':
    main()
