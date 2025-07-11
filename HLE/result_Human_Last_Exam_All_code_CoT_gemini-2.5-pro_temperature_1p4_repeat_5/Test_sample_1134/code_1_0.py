import numpy as np

def ket_to_str(ket, precision=3):
    """Converts a ket vector to a readable string format."""
    ket = np.round(ket, precision)
    if np.allclose(ket, [1, 0]): return "|0>"
    if np.allclose(ket, [0, 1]): return "|1>"
    if np.allclose(ket / (1/np.sqrt(2)), [1, 1]): return "|+>"
    if np.allclose(ket / (1/np.sqrt(2)), [1, -1]): return "|->"
    if np.allclose(ket / (1/np.sqrt(2)), [1, 1j]): return "|i>"
    if np.allclose(ket / (1/np.sqrt(2)), [1, -1j]): return "|-i>"
    # Normalize for phase check
    if np.isclose(np.linalg.norm(ket), 1.0):
        if np.allclose(np.abs(ket), [0,1]): return ket[1]/1j*np.exp(1j*np.pi/2) and "|1>" # phase on |1>
    return f"[{ket[0]:.2f}, {ket[1]:.2f}]"

# --- State Definitions ---
ket0 = np.array([1, 0], dtype=complex)
ket1 = np.array([0, 1], dtype=complex)
s_minus = (ket0 - ket1) / np.sqrt(2)
s_i = (ket0 + 1j * ket1) / np.sqrt(2)
s_neg_i = (ket0 - 1j * ket1) / np.sqrt(2)

initial_states = {
    "|0>": ket0,
    "|1>": ket1,
    "|->": s_minus,
    "|i>": s_i,
    "|-i>": s_neg_i
}

# --- Gate Definition ---
# C. S dagger gate
S_dag = np.array([[1, 0], [0, -1j]], dtype=complex)
gate_name = "S†"
gate_matrix = S_dag

print(f"Analyzing the {gate_name} gate:")
failure_found = False

for name, state_ket in initial_states.items():
    # Apply the gate
    final_state_ket = gate_matrix @ state_ket
    
    # Check if the result is a death state
    if np.allclose(final_state_ket, s_i) or np.allclose(final_state_ket, s_neg_i):
        status = "UNSAFE"
        failure_found = True
    else:
        status = "Safe"
        
    print(f"  {gate_name} {name} = {ket_to_str(final_state_ket)}  ...  ({status})")

if failure_found:
    print("\nThe action is not perfectly safe. The critical equation showing the failure is:")
    # The matrix representation of the equation S_dag * |-> = |i>
    # [1, 0]   [1/sqrt(2)]   [ 1/sqrt(2)]
    # [0,-i] * [-1/sqrt(2)] = [ i/sqrt(2)]
    print(f"{gate_name} |-> = |i>")
    print("Which corresponds to the calculation:")
    print("[[1, 0]  [ 0.707]]   [ 0.707]")
    print("[[0,-i]] * [-0.707]] = [ i*0.707]")
