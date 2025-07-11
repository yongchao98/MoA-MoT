import cmath

def solve_problem():
    """
    Calculates the rank of H^2_c(Y, Q) for the given algebraic variety Y.
    """
    # The orientation-preserving icosahedral group is isomorphic to A_5.
    # We consider its standard 3D representation, which is a subgroup of SO(3).
    # The problem is to find the number of "junior" conjugacy classes.
    # A class is junior if the "age" of its elements is 1.
    # age(g) = sum(theta_j) where eigenvalues are lambda_j = exp(2*pi*i*theta_j)
    # with 0 <= theta_j < 1.

    # The character table of A_5 gives the trace of elements in each class for the 3D representation.
    # Let phi = (1 + sqrt(5)) / 2.
    phi = (1 + 5**0.5) / 2
    conjugacy_classes = [
        {'name': 'Identity', 'order': 1, 'char': 3},
        {'name': 'Rotation by pi', 'order': 2, 'char': -1},
        {'name': 'Rotation by 2pi/3', 'order': 3, 'char': 0},
        {'name': 'Rotation by 2pi/5', 'order': 5, 'char': phi},
        {'name': 'Rotation by 4pi/5', 'order': 5, 'char': 1 - phi},
    ]

    def get_phases_from_eigenvalues(eigenvalues):
        """Calculates phases theta in [0, 1) for eigenvalues exp(2*pi*i*theta)."""
        phases = []
        for lmbda in eigenvalues:
            # cmath.phase is in (-pi, pi], we want a result in [0, 1) after division.
            phase = cmath.phase(lmbda) / (2 * cmath.pi)
            if phase < 0:
                phase += 1
            # Handle potential floating point inaccuracies close to 1
            if abs(phase - 1.0) < 1e-9:
                phase = 0.0
            phases.append(phase)
        # Sort for consistent output.
        phases.sort()
        return phases

    junior_classes_count = 0

    print("The rank of H^2_c(Y, Q) is equal to b_4(Y), which for a crepant resolution")
    print("of C^3/G is equal to b_2(Y). By the McKay correspondence, b_2(Y) is the")
    print("number of 'junior' conjugacy classes (age = 1) of the group G = Icosahedral group (A_5).\n")
    print("Calculating the age for each conjugacy class:")
    print("---------------------------------------------")

    # --- Class 1: Identity ---
    c = conjugacy_classes[0]
    print(f"Analyzing class '{c['name']}' (order {c['order']}, character {c['char']}):")
    # For an element g in SO(3), eigenvalues are 1, e^{i*alpha}, e^{-i*alpha}
    # Identity has alpha = 0. Eigenvalues are (1, 1, 1).
    eigenvalues = [1, 1, 1]
    phases = get_phases_from_eigenvalues(eigenvalues)
    age = sum(phases)
    print(f"  Eigenvalues are (1, 1, 1).")
    print(f"  Phases are ({phases[0]:.3f}, {phases[1]:.3f}, {phases[2]:.3f}).")
    print(f"  Age = {phases[0]:.3f} + {phases[1]:.3f} + {phases[2]:.3f} = {age:.2f}")
    if abs(age - 1.0) < 1e-9:
        junior_classes_count += 1
        print("  This class is junior.\n")
    else:
        print("  This class is NOT junior.\n")

    # --- Class 2: Order 2 ---
    c = conjugacy_classes[1]
    print(f"Analyzing class '{c['name']}' (order {c['order']}, character {c['char']}):")
    # Rotation by pi. Eigenvalues are (1, -1, -1). Sum is -1, matches character.
    eigenvalues = [1, -1, -1]
    phases = get_phases_from_eigenvalues(eigenvalues)
    age = sum(phases)
    print(f"  Eigenvalues are (1, -1, -1).")
    print(f"  Phases are ({phases[0]:.3f}, {phases[1]:.3f}, {phases[2]:.3f}).")
    print(f"  Age = {phases[0]:.3f} + {phases[1]:.3f} + {phases[2]:.3f} = {age:.2f}")
    if abs(age - 1.0) < 1e-9:
        junior_classes_count += 1
        print("  This class is junior.\n")
    else:
        print("  This class is NOT junior.\n")

    # --- Class 3: Order 3 ---
    c = conjugacy_classes[2]
    print(f"Analyzing class '{c['name']}' (order {c['order']}, character {c['char']}):")
    # Rotation by 2pi/3. Eigenvalues are (1, e^{2pi*i/3}, e^{-2pi*i/3})
    # Sum is 1 + w + w^2 = 0, which matches the character.
    w = cmath.exp(2j * cmath.pi / 3)
    eigenvalues = [1, w, w.conjugate()]
    phases = get_phases_from_eigenvalues(eigenvalues)
    age = sum(phases)
    print(f"  Eigenvalues are (1.00, {w.real:.2f}{w.imag:+.2f}j, {w.conjugate().real:.2f}{w.conjugate().imag:+.2f}j).")
    print(f"  Phases are ({phases[0]:.3f}, {phases[1]:.3f}, {phases[2]:.3f}).")
    print(f"  Age = {phases[0]:.3f} + {phases[1]:.3f} + {phases[2]:.3f} = {age:.2f}")
    if abs(age - 1.0) < 1e-9:
        junior_classes_count += 1
        print("  This class is junior.\n")
    else:
        print("  This class is NOT junior.\n")

    # --- Class 4: Order 5 (char = phi) ---
    c = conjugacy_classes[3]
    print(f"Analyzing class '{c['name']}' (order {c['order']}, character {c['char']:.3f}):")
    # Rotation by 2pi/5. Eigenvalues are (1, e^{2pi*i/5}, e^{-2pi*i/5})
    # Sum is 1 + 2*cos(2pi/5) = phi.
    v = cmath.exp(2j * cmath.pi / 5)
    eigenvalues = [1, v, v.conjugate()]
    phases = get_phases_from_eigenvalues(eigenvalues)
    age = sum(phases)
    print(f"  Eigenvalues are (1.00, {v.real:.2f}{v.imag:+.2f}j, {v.conjugate().real:.2f}{v.conjugate().imag:+.2f}j).")
    print(f"  Phases are ({phases[0]:.3f}, {phases[1]:.3f}, {phases[2]:.3f}).")
    print(f"  Age = {phases[0]:.3f} + {phases[1]:.3f} + {phases[2]:.3f} = {age:.2f}")
    if abs(age - 1.0) < 1e-9:
        junior_classes_count += 1
        print("  This class is junior.\n")
    else:
        print("  This class is NOT junior.\n")

    # --- Class 5: Order 5 (char = 1-phi) ---
    c = conjugacy_classes[4]
    print(f"Analyzing class '{c['name']}' (order {c['order']}, character {c['char']:.3f}):")
    # Rotation by 4pi/5. Eigenvalues are (1, e^{4pi*i/5}, e^{-4pi*i/5})
    # Sum is 1 + 2*cos(4pi/5) = 1-phi.
    u = cmath.exp(4j * cmath.pi / 5)
    eigenvalues = [1, u, u.conjugate()]
    phases = get_phases_from_eigenvalues(eigenvalues)
    age = sum(phases)
    print(f"  Eigenvalues are (1.00, {u.real:.2f}{u.imag:+.2f}j, {u.conjugate().real:.2f}{u.conjugate().imag:+.2f}j).")
    print(f"  Phases are ({phases[0]:.3f}, {phases[1]:.3f}, {phases[2]:.3f}).")
    print(f"  Age = {phases[0]:.3f} + {phases[1]:.3f} + {phases[2]:.3f} = {age:.2f}")
    if abs(age - 1.0) < 1e-9:
        junior_classes_count += 1
        print("  This class is junior.\n")
    else:
        print("  This class is NOT junior.\n")

    print("---------------------------------------------")
    print(f"Total number of junior classes found: {junior_classes_count}")
    print(f"The rank of H^2_c(Y, Q) is therefore {junior_classes_count}.")

solve_problem()