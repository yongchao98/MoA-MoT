import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.interpolate import interp1d

def potential(x, alpha):
    """Defines the sextic anharmonic oscillator potential."""
    return 0.5 * (alpha**2 - 7) * x**2 - alpha * x**4 + 0.5 * x**6

def schrodinger_ode(x, y, E, alpha):
    """Defines the ODE system for the time-independent Schrodinger equation."""
    # y = [psi, dpsi/dx]
    psi, dpsi = y
    d2psi = 2 * (potential(x, alpha) - E) * psi
    return [dpsi, d2psi]

def shoot(E, alpha, x_max=10.0):
    """
    Solves the ODE for a given energy E and returns the value of psi at x_max.
    We solve for even parity states: psi(0)=1, psi'(0)=0.
    """
    sol = solve_ivp(
        schrodinger_ode,
        [0, x_max],
        [1.0, 0.0],  # Initial conditions for even parity
        args=(E, alpha),
        dense_output=True,
        t_eval=[x_max]
    )
    return sol.y[0, -1]

def find_eigenvalue(alpha, n_even, x_max=10.0, E_step=0.1, max_steps=2000):
    """
    Finds the n-th even eigenvalue for a given alpha using the shooting method.
    n_even = 0 for ground state (E0), 1 for second excited state (E2).
    """
    # Calculate the minimum of the potential for a good energy search start
    # Extrema are at x=0 and x_ext^2 = (2*alpha + sqrt(alpha^2 + 21)) / 3
    x_ext_sq = (2 * alpha + np.sqrt(alpha**2 + 21)) / 3.0
    v_at_x_ext = potential(np.sqrt(x_ext_sq), alpha)
    v_at_0 = potential(0, alpha)  # This is always 0
    E_start = min(v_at_x_ext, v_at_0) - 1.0  # Start slightly below V_min

    roots_found = 0
    E_low = E_start
    val_low = shoot(E_low, alpha, x_max)

    for i in range(max_steps):
        E_high = E_low + E_step
        val_high = shoot(E_high, alpha, x_max)
        
        # Check for a sign change, indicating an eigenvalue in the bracket
        if np.sign(val_low) != np.sign(val_high):
            if roots_found == n_even:
                # Bracket found for the desired eigenvalue, find it precisely
                eigenvalue = brentq(shoot, E_low, E_high, args=(alpha, x_max))
                return eigenvalue
            roots_found += 1
        
        E_low = E_high
        val_low = val_high
    
    raise RuntimeError(f"Could not find the eigenvalue for n_even={n_even}, alpha={alpha}")

def get_psi2_interpolated(alpha, x_max=10.0, num_points=500):
    """
    Calculates the normalized second excited state wavefunction psi_2(x)
    for a given alpha and returns it as an interpolated function.
    """
    # Find the second even eigenvalue, E_2 (n_even=1)
    E2 = find_eigenvalue(alpha, n_even=1, x_max=x_max)

    # Solve the ODE again with the correct eigenvalue E2 to get the wavefunction
    sol = solve_ivp(
        schrodinger_ode,
        [0, x_max],
        [1.0, 0.0],
        args=(E2, alpha),
        dense_output=True,
        t_eval=np.linspace(0, x_max, num_points)
    )
    
    # Normalize the wavefunction (integrate psi^2 from -x_max to x_max)
    # Since psi is even, integral is 2 * integral from 0 to x_max.
    psi = sol.y[0]
    x_dense = sol.t
    norm_sq = np.trapz(psi**2, x_dense)
    norm = np.sqrt(2 * norm_sq)
    psi_normalized = psi / norm
    
    # Return an interpolated function (valid for x >= 0)
    return interp1d(x_dense, psi_normalized, bounds_error=False, fill_value=0)

def target_function(alpha):
    """
    The function g(alpha) = psi_2(alpha; x=alpha) we want to find the root of.
    """
    try:
        # Get the interpolated function for psi_2(x) for this alpha
        psi2_func = get_psi2_interpolated(alpha)
        # Evaluate psi_2 at x=alpha
        value = psi2_func(alpha)
        return value
    except RuntimeError as e:
        # If eigenvalue search fails, return a large value to guide root finder away
        print(f"Warning: {e}")
        return 1e6

if __name__ == "__main__":
    # Bracket for alpha where the root is expected to lie.
    alpha_low = 2.5
    alpha_high = 3.5

    print(f"Searching for alpha_0 in the interval [{alpha_low}, {alpha_high}]...")

    try:
        # Find the root alpha_0 where target_function(alpha) = 0
        alpha_0 = brentq(target_function, alpha_low, alpha_high, xtol=1e-7)
        
        # The final equation is F(alpha_0) = 0.
        # As requested, here is the number in that equation.
        print("Final value alpha_0 found:")
        print(alpha_0)
        
    except ValueError:
        print(f"Failed to find a root. The function values at the endpoints "
              f"f({alpha_low})={target_function(alpha_low):.4f} and "
              f"f({alpha_high})={target_function(alpha_high):.4f} "
              f"may not have opposite signs.")
    except Exception as e:
        print(f"An error occurred: {e}")
