import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import warnings

# Suppress warnings that may arise from the solver for non-eigenvalue energies
warnings.filterwarnings("ignore", category=UserWarning)

# --- Step 1: Define the physical system ---

def potential(x, alpha):
    """The sextic anharmonic oscillator potential V(x)."""
    return -3.5 * x**2 + 0.5 * alpha**2 * x**2 - alpha * x**4 + 0.5 * x**6

def schrodinger_ode(x, y, E, alpha):
    """The time-independent Schrödinger equation as a system of 1st order ODEs."""
    psi, phi = y  # y = [psi, psi']
    return [phi, 2 * (potential(x, alpha) - E) * psi]

# --- Step 2: Implement the numerical solver (shooting method) ---

def objective_function(E, alpha, n, x_max):
    """
    Objective function for the root finder.
    It returns the value of the wavefunction at a large distance x_max.
    Eigenvalues are the energies E for which this function is zero.
    """
    is_even = n % 2 == 0
    y0 = [1.0, 0.0] if is_even else [0.0, 1.0] # Initial conditions for even/odd states
    
    sol = solve_ivp(
        schrodinger_ode,
        [0, x_max],
        y0,
        args=(E, alpha),
        t_eval=[x_max] # We only need the final value for the root finder
    )
    return sol.y[0, -1]

def find_eigenstate(n, alpha, x_max_base=10.0, E_max_search=200.0, n_scan=500):
    """
    Finds the n-th eigenenergy and eigenfunction for a given alpha.
    It works by finding the roots of the objective_function.
    """
    is_even = n % 2 == 0
    if not is_even:
        raise NotImplementedError("This solver is configured only for even states (n=0, 2, ...).")

    # The n-th even state corresponds to the (n/2)-th root of the objective function.
    target_root_index = n // 2
    
    # Adjust integration range based on alpha
    x_max = max(x_max_base, alpha + 3.0)

    # Estimate the potential minimum to set a search range for energy
    search_x = np.linspace(0, x_max, 100)
    try:
        v_min = np.min(potential(search_x, alpha))
    except (ValueError, TypeError):
        return None, None, None

    # Scan for energy brackets that contain eigenvalues
    E_scan = np.linspace(v_min, E_max_search, n_scan)
    psi_ends = np.array([objective_function(E, alpha, n, x_max) for E in E_scan])
    
    # Find indices before a sign change, indicating a bracketed root
    sign_changes = np.where(np.diff(np.sign(psi_ends)))[0]

    if len(sign_changes) <= target_root_index:
        # Failed to find enough eigenvalues in the search range
        return None, None, None

    # Get the energy bracket for the desired state
    bracket_index = sign_changes[target_root_index]
    E_low = E_scan[bracket_index]
    E_high = E_scan[bracket_index + 1]

    try:
        # Find the precise eigenvalue using the bracket
        eigen_E = brentq(objective_function, E_low, E_high, args=(alpha, n, x_max))
        
        # Now, solve the ODE with the correct energy to get the full wavefunction
        sol = solve_ivp(
            schrodinger_ode,
            [0, x_max],
            [1.0, 0.0], # Initial conditions for an even state
            args=(eigen_E, alpha),
            dense_output=True,
            t_eval=np.linspace(0, x_max, 1001)
        )
        return eigen_E, sol.t, sol.y[0]
    except (ValueError, RuntimeError):
        # brentq can fail if the bracket is not valid
        return None, None, None

# --- Step 3: Define the functions whose roots we need to find ---

def find_alpha_for_E2_zero(alpha):
    """Function wrapper for root finding. Returns E_2(alpha)."""
    if alpha <= 0: return 1e9
    result = find_eigenstate(n=2, alpha=alpha)
    if result is None or result[0] is None:
        return 1e9  # Return a large value on failure
    return result[0] # Return E_2

def find_alpha_for_psi2_zero(alpha):
    """Function wrapper for root finding. Returns psi_2(alpha; alpha)."""
    if alpha <= 0: return 1e9
    E2, x, psi2 = find_eigenstate(n=2, alpha=alpha)
    if E2 is None:
        return 1e9 # Return a large value on failure

    # Interpolate to find psi_2 at x=alpha
    if alpha > x[-1]: # Check if alpha is within interpolation range
         return 1e9
    psi2_interp = interp1d(x, psi2, kind='cubic', bounds_error=False, fill_value=1e9)
    return float(psi2_interp(alpha))

# --- Step 4: Search for the roots and find the largest value ---

def solve_problem():
    """
    Main function to find all roots for F(alpha)=0 and determine the largest.
    """
    print("Searching for alpha_0 where F(alpha_0) = 0...")
    print("This requires solving the Schrödinger equation numerically and may take a moment.")
    
    all_roots = []
    
    # Define the scan range for alpha
    alpha_scan = np.linspace(0.1, 7.0, 150)

    # --- Find roots for E_2(alpha) = 0 ---
    print("\nSearching for roots of E_2(alpha) = 0...")
    f_vals_E2 = [find_alpha_for_E2_zero(a) for a in alpha_scan]
    sign_changes_E2 = np.where(np.diff(np.sign(f_vals_E2)))[0]
    
    for idx in sign_changes_E2:
        a_low, a_high = alpha_scan[idx], alpha_scan[idx+1]
        try:
            root = brentq(find_alpha_for_E2_zero, a_low, a_high, xtol=1e-9)
            all_roots.append({'value': root, 'condition': 'E_2(alpha) = 0'})
            print(f"Found root alpha = {root:.8f} from E_2(alpha) = 0")
        except (ValueError, RuntimeError):
            print(f"Could not refine root in bracket [{a_low:.2f}, {a_high:.2f}] for E_2")

    # --- Find roots for psi_2(alpha, alpha) = 0 ---
    print("\nSearching for roots of psi_2(alpha; alpha) = 0...")
    f_vals_psi2 = [find_alpha_for_psi2_zero(a) for a in alpha_scan]
    sign_changes_psi2 = np.where(np.diff(np.sign(f_vals_psi2)))[0]

    for idx in sign_changes_psi2:
        a_low, a_high = alpha_scan[idx], alpha_scan[idx+1]
        try:
            root = brentq(find_alpha_for_psi2_zero, a_low, a_high, xtol=1e-9)
            all_roots.append({'value': root, 'condition': 'psi_2(alpha; alpha) = 0'})
            print(f"Found root alpha = {root:.8f} from psi_2(alpha; alpha) = 0")
        except (ValueError, RuntimeError):
            print(f"Could not refine root in bracket [{a_low:.2f}, {a_high:.2f}] for psi_2")
            
    if not all_roots:
        print("\nNo roots found in the specified range.")
        return None

    # --- Find the largest root ---
    largest_root_info = max(all_roots, key=lambda r: r['value'])
    alpha_0 = largest_root_info['value']
    condition = largest_root_info['condition']

    print(f"\n--------------------------------------------------")
    print(f"The largest value alpha_0 such that F(alpha_0) = 0 is: {alpha_0:.8f}")
    print(f"This root satisfies the condition: {condition}")
    
    # Fulfilling the request to "output each number in the final equation"
    if 'E_2' in condition:
        # Equation is E_2(alpha_0) = 0
        print(f"The numbers in the final equation E_2({alpha_0:.8f}) = 0 are 2, {alpha_0:.8f}, and 0.")
    elif 'psi_2' in condition:
        # Equation is psi_2(alpha_0; alpha_0) = 0
        print(f"The numbers in the final equation psi_2({alpha_0:.8f}; {alpha_0:.8f}) = 0 are 2, {alpha_0:.8f}, and 0.")
        
    return alpha_0

# Execute the solver and store the final answer
alpha_final = solve_problem()
if alpha_final is not None:
    print(f"\nFinal Answer: {alpha_final}")
