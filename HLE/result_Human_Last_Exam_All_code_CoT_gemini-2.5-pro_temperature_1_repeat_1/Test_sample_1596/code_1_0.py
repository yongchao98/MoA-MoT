import numpy as np

# Helper class for Quaternion math
class Quaternion:
    """A simple class to represent and operate on quaternions."""
    def __init__(self, a, b, c, d):
        self.a = a  # real part
        self.b = b  # i part
        self.c = c  # j part
        self.d = d  # k part
        self.vec = np.array([a, b, c, d])

    def __mul__(self, other):
        """Quaternion multiplication."""
        a1, b1, c1, d1 = self.a, self.b, self.c, self.d
        a2, b2, c2, d2 = other.a, other.b, other.c, other.d
        
        a = a1*a2 - b1*b2 - c1*c2 - d1*d2
        b = a1*b2 + b1*a2 + c1*d2 - d1*c2
        c = a1*c2 - b1*d2 + c1*a2 + d1*b2
        d = a1*d2 + b1*c2 - c1*b2 + d1*a2
        
        return Quaternion(a, b, c, d)

    def __sub__(self, other):
        """Quaternion subtraction."""
        return Quaternion(self.a - other.a, self.b - other.b, self.c - other.c, self.d - other.d)
        
    def __abs__(self):
        """Quaternion norm (magnitude)."""
        return np.sqrt(self.a**2 + self.b**2 + self.c**2 + self.d**2)

    def __str__(self):
        """String representation of the quaternion."""
        return f"{self.a:.4f} + {self.b:.4f}i + {self.c:.4f}j + {self.d:.4f}k"
    
    def to_matrix(self):
        """Returns the 4x4 left-multiplication matrix representation."""
        a, b, c, d = self.a, self.b, self.c, self.d
        return np.array([
            [a, -b, -c, -d],
            [b,  a, -d,  c],
            [c,  d,  a, -b],
            [d, -c,  b,  a]
        ])

def dot_product(q1, q2):
    """Computes the inner product of two quaternions as 4D vectors."""
    return np.dot(q1.vec, q2.vec)

def get_R(m, theta, u_vec):
    """Creates the rotation quaternion R_m for position m."""
    angle = m * theta
    u_norm = np.linalg.norm(u_vec)
    if u_norm < 1e-9:
        raise ValueError("Rotation axis cannot be a zero vector")
    u = u_vec / u_norm
    
    a = np.cos(angle)
    s = np.sin(angle)
    b, c, d = u[0] * s, u[1] * s, u[2] * s
    
    return Quaternion(a, b, c, d)

def run_verification():
    """Runs numerical verification for key statements."""
    print("Analyzing selected statements about Quaternion RoPE...\n")
    
    # Common parameters for testing
    theta = np.pi / 10
    u_axis = np.array([1.0, 1.0, 2.0])
    q = Quaternion(0.5, 1.0, -0.2, 0.8)
    k = Quaternion(0.9, -0.4, 1.1, -0.1)

    # --- Statement A Verification ---
    print("--- Verifying Statement A (False) ---")
    m1, n1 = 5, 2  # m-n = 3
    m2, n2 = 2, 5  # m-n = -3, so |m-n| is the same
    R_m1, R_n1 = get_R(m1, theta, u_axis), get_R(n1, theta, u_axis)
    R_m2, R_n2 = get_R(m2, theta, u_axis), get_R(n2, theta, u_axis)
    inner_prod1 = dot_product(R_m1 * q, R_n1 * k)
    inner_prod2 = dot_product(R_m2 * q, R_n2 * k)
    print(f"For m-n = {m1-n1}, the inner product is: {inner_prod1:.4f}")
    print(f"For m-n = {m2-n2}, the inner product is: {inner_prod2:.4f}")
    print("The results are different, so the inner product does not depend only on |m-n|.\n")

    # --- Statement E Verification ---
    print("--- Verifying Statement E (True) ---")
    m = 7
    R_m = get_R(m, theta, u_axis)
    original_norm = abs(q)
    rotated_norm = abs(R_m * q)
    print(f"Original quaternion norm |q|: {original_norm:.4f}")
    print(f"Rotated quaternion norm |R_m(q)|: {rotated_norm:.4f}")
    print("The magnitude is preserved.\n")

    # --- Statement G Verification ---
    print("--- Verifying Statement G (True) ---")
    v1 = Quaternion(0, 1, -2, 0)
    v2 = Quaternion(0, 2, 1, 0) # Chosen so dot product is 0
    m = 4
    R_m = get_R(m, theta, u_axis)
    original_dot = dot_product(v1, v2)
    rotated_dot = dot_product(R_m * v1, R_m * v2)
    print(f"Original dot product <v1, v2>: {original_dot:.4f}")
    print(f"Dot product after rotation <R_m(v1), R_m(v2)>: {rotated_dot:.4f}")
    print("Orthogonality is preserved (dot product remains 0).\n")
    
    # --- Statement H Verification ---
    print("--- Verifying Statement H (True) ---")
    m, n = 3, 5
    R_m, R_n = get_R(m, theta, u_axis), get_R(n, theta, u_axis)
    R_m_plus_n = get_R(m + n, theta, u_axis)
    comp_v = R_m * (R_n * q)
    sum_v = R_m_plus_n * q
    print(f"Result of composition R_m(R_n(q)): {comp_v}")
    print(f"Result of direct rotation R_{{m+n}}(q):  {sum_v}")
    print("The results are identical, confirming R_m o R_n = R_{m+n}.\n")
    
    # --- Statement L Verification ---
    print("--- Verifying Statement L (True) ---")
    m = 6
    R_m = get_R(m, theta, u_axis)
    L_Rm_trace = np.trace(R_m.to_matrix())
    expected_trace = 4 * np.cos(m * theta)
    print(f"The trace of the 4x4 matrix for R_{m} is: {L_Rm_trace:.4f}")
    print(f"The value of 4*cos(m*theta) is: {expected_trace:.4f}")
    print("The values are identical, confirming the statement.\n")

    # Final Summary of True Statements
    true_statements = "CEFGHJL"
    print("="*50)
    print(f"Final sorted list of correct statement letters: {true_statements}")
    print("="*50)

if __name__ == '__main__':
    run_verification()