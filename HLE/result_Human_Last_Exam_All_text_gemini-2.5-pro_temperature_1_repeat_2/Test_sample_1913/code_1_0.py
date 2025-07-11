# E = EllipticCurve("161083.a1") or EllipticCurve([0, -1, 1, -10, -20])
# K.<z3> = CyclotomicField(3)
# G = DirichletGroup(7, K)
# chi = [c for c in G if c.order() == 3][0]
# L_twist = E.lseries().twist(chi)
# a = L_twist.derivative(1)