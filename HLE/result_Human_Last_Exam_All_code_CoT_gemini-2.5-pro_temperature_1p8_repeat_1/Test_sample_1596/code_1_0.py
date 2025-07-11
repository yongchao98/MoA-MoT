import numpy as np

class Quaternion:
    """A class to represent and operate on quaternions."""
    def __init__(self, w, x=None, y=None, z=None):
        if isinstance(w, (list, tuple, np.ndarray)):
            self.coeffs = np.array(w, dtype=float)
        else:
            self.coeffs = np.array([w, x, y, z], dtype=float)

    @property
    def w(self): return self.coeffs[0]
    @property
    def x(self): return self.coeffs[1]
    @property
    def y(self): return self.coeffs[2]
    @property
    def z(self): return self.coeffs[3]
    @property
    def re(self): return self.w
    @property
    def vec(self): return self.coeffs[1:]

    def __repr__(self):
        return f"({self.w:.3f} + {self.x:.3f}i + {self.y:.3f}j + {self.z:.3f}k)"

    def __mul__(self, other):
        if not isinstance(other, Quaternion):
            if isinstance(other, (int, float)):
                return Quaternion(self.coeffs * other)
            return NotImplemented
        
        a1, b1, c1, d1 = self.coeffs
        a2, b2, c2, d2 = other.coeffs

        w = a1*a2 - b1*b2 - c1*c2 - d1*d2
        x = a1*b2 + b1*a2 + c1*d2 - d1*c2
        y = a1*c2 - b1*d2 + c1*a2 + d1*b2
        z = a1*d2 + b1*c2 - c1*b2 + d1*a2
        return Quaternion(w, x, y, z)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        return Quaternion(self.coeffs + other.coeffs)

    def __sub__(self, other):
        return Quaternion(self.coeffs - other.coeffs)

    def norm(self):
        return np.linalg.norm(self.coeffs)

    def __eq__(self, other):
        return np.allclose(self.coeffs, other.coeffs)

def inner_product(p: Quaternion, q: Quaternion) -> float:
    """Computes the inner product of two quaternions."""
    return np.dot(p.coeffs, q.coeffs)

def get_Rm(m: int, theta: float, u_vec: np.ndarray) -> Quaternion:
    """Generates the RoPE rotation quaternion R_m."""
    u_norm = np.linalg.norm(u_vec)
    if not np.isclose(u_norm, 1.0):
       u_vec = u_vec / u_norm
    
    angle = m * theta
    w = np.cos(angle)
    s = np.sin(angle)
    return Quaternion(w, *(s * u_vec))

def solve():
    """
    Analyzes each statement and returns the sorted string of true statements.
    """
    # Setup for numerical checks
    theta = np.pi / 8
    u_vec = np.array([3, 4, 12]) / 13.0
    u_quat = Quaternion(0, *u_vec)
    q = Quaternion(1, 2, -1, 0.5)
    k = Quaternion(-0.5, 1, 3, -2)
    v = Quaternion(2, -1, 0.5, 1.5)
    unit_v = Quaternion(v.coeffs / v.norm())
    m, n, p = 5, 3, 2
    alpha = 2.5
    
    true_statements = []

    # Statement A: ⟨R_m(q), R_n(k)⟩ depends only on |m-n|
    # Analysis: False. It depends on (m-n). ⟨q,R_{m-n}*k⟩ vs ⟨q,R_{n-m}*k⟩ are different.
    m1, n1 = 5, 3 # diff = 2
    m2, n2 = 3, 5 # diff = -2
    ip1 = inner_product(get_Rm(m1, theta, u_vec) * q, get_Rm(n1, theta, u_vec) * k)
    ip2 = inner_product(get_Rm(m2, theta, u_vec) * q, get_Rm(n2, theta, u_vec) * k)
    if not np.allclose(ip1, ip2):
        pass # A is False
    else:
        true_statements.append('A')

    # Statement B: R_m(q)R_n(k) = R_{m+p}(q)R_{n+p}(k) for any p
    # Analysis: False. Quaternion multiplication is not commutative with rotation like this.
    lhs_B = (get_Rm(m, theta, u_vec) * q) * (get_Rm(n, theta, u_vec) * k)
    rhs_B = (get_Rm(m + p, theta, u_vec) * q) * (get_Rm(n + p, theta, u_vec) * k)
    if not (lhs_B == rhs_B):
        pass # B is False
    else:
        true_statements.append('B')
        
    # Statement C: The rotation axis (u_1,u_2,u_3) must be fixed for all positions
    # Analysis: True. To get the relative property ⟨R_m(q),R_n(k)⟩ = ⟨R_{m-n}(q),k⟩, the axis must be constant.
    true_statements.append('C')
    
    # Statement D: Quaternion RoPE can encode 3D relative positions with a single rotation
    # Analysis: False. The formulation uses a scalar position `m`, so it encodes 1D positions.
    pass # D is False

    # Statement E: The magnitude |R_m(v)| equals |v| for all m
    # Analysis: True. R_m is a unit quaternion, and |R_m*v| = |R_m|*|v| = 1*|v|.
    if np.allclose((get_Rm(m, theta, u_vec) * v).norm(), v.norm()):
        true_statements.append('E')
    
    # Statement F: R_m(αv) = αR_m(v) for scalar α
    # Analysis: True. Scalar multiplication is commutative in quaternion algebra.
    lhs_F = get_Rm(m, theta, u_vec) * (alpha * v)
    rhs_F = alpha * (get_Rm(m, theta, u_vec) * v)
    if lhs_F == rhs_F:
        true_statements.append('F')
        
    # Statement G: The quaternion rotation preserves orthogonality between vectors
    # Analysis: True. Left multiplication by a unit quaternion is an isometry, so ⟨Rp,Rq⟩ = ⟨p,q⟩.
    p_orth = Quaternion(1, 1, 0, 0)
    q_orth = Quaternion(0, 0, 1, 1) # inner_product is 0
    if np.allclose(inner_product(p_orth, q_orth), inner_product(get_Rm(m, theta, u_vec)*p_orth, get_Rm(m, theta, u_vec)*q_orth)):
        true_statements.append('G')

    # Statement H: The composition R_m ∘ R_n equals R_{m+n}
    # Analysis: True. This is De Moivre's formula for quaternions with the same rotation axis.
    if get_Rm(m, theta, u_vec) * get_Rm(n, theta, u_vec) == get_Rm(m + n, theta, u_vec):
        true_statements.append('H')

    # Statement J: (R_m ∘ R_n)(v) - R_n(R_m(v)) is always purely imaginary
    # Analysis: True. R_m and R_n commute since they share an axis, so the result is the zero quaternion, which has a zero real part.
    commutator_R = get_Rm(m, theta, u_vec) * get_Rm(n, theta, u_vec) - get_Rm(n, theta, u_vec) * get_Rm(m, theta, u_vec)
    if np.allclose((commutator_R * v).re, 0):
        true_statements.append('J')

    # Statement K: R_m(R_n(v)) preserves the real component of v regardless of m,n
    # Analysis: False. Re(R_{m+n}v) = v_0*cos - sin*(u.v_vec), which is not generally v_0.
    v_after = get_Rm(m+n, theta, u_vec) * v
    if not np.allclose(v.re, v_after.re):
        pass # K is False
    else:
        #This might be true for specific values, need to check with another
        v2 = Quaternion(1,2,3,4)
        v_after_2 = get_Rm(m+n, theta, u_vec) * v2
        if not np.allclose(v2.re, v_after_2.re):
             pass
        else:
             true_statements.append('K')
             
    # Statement L: The trace of the 4×4 matrix representation of R_m equals 4cos(mθ)
    # Analysis: True. The trace of the left-multiplication matrix for q is 4*Re(q). Re(R_m) is cos(mθ).
    if np.allclose(4 * get_Rm(m, theta, u_vec).re, 4 * np.cos(m * theta)):
        true_statements.append('L')

    # Statement M: The commutator [R_m, R_n] = 0 if and only if the rotation axis is parallel to (1,1,1)
    # Analysis: False. The commutator is always zero for this formulation, regardless of the axis. The "only if" makes it false.
    pass # M is False
    
    # Statement N: For any unit quaternion v, lim_{m→∞} |R_m(v)| converges to 1/√2
    # Analysis: False. The norm is always preserved, so |R_m(v)| = |v| = 1. The limit is 1.
    pass # N is False
    
    print("".join(sorted(true_statements)))

solve()
<<<CEFGHJL>>>